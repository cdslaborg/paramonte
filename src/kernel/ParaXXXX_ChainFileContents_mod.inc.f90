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

!>  \brief This module contains the classs and procedures for chain IO and manipulation.
!>  \author Amir Shahmoradi

#if defined PARADRAM
#define ParaXXXX ParaDRAM
#elif defined PARADISE
#define ParaXXXX ParaDISE
#elif defined PARANEST
#define ParaXXXX ParaNest
#else
#error "Unrecognized sampler in ParaXXXX_ChainFileContents_mod.inc.f90"
#endif

    use, intrinsic :: iso_fortran_env, only: output_unit
    use JaggedArray_mod, only: CharVec_type
    use Decoration_mod, only: INDENT
    use Constants_mod, only: IK, RK
    use Constants_mod, only: PMSM
    use Err_mod, only: Err_type
    use Err_mod, only: warn
    implicit none

    character(*), parameter :: MODULE_NAME = "@"//PMSM%ParaXXXX//"@ParaMonteChainFileContents_mod"

#if defined PARADRAM || defined PARADISE

    type :: Count_type
        integer(IK) :: verbose = 0_IK   ! the number of points (weight=1) in the chain.
        integer(IK) :: compact = 0_IK   ! the number of unique (weighted) points in the chain.
        integer(IK) :: target = 0_IK    ! the size of the allocations for the Chain components.
    end type Count_type

    character(*), parameter :: COL_HEADER_DEFAULT(*) =  [ "ProcessID            " &
                                                        , "DelayedRejectionStage" &
                                                        , "MeanAcceptanceRate   " &
                                                        , "AdaptationMeasure    " &
                                                        , "BurninLocation       " &
                                                        , "SampleWeight         " &
                                                        , "SampleLogFunc        " &
                                                        ]

#elif defined PARANEST

    type :: Count_type
        real)RK)    :: verbose = 0._RK  ! the number of points (weight=1) in the chain.
        integer(IK) :: compact = 0_IK   ! the number of unique (weighted) points in the chain.
        integer(IK) :: target = 0_IK    ! the size of the allocations for the Chain components.
    end type Count_type

    character(*), parameter :: COL_HEADER_DEFAULT(*) =  [ "ProcessID            " &
                                                        , "MeanAcceptanceRate   " &
                                                        , "RemainingPriorMass   " &
                                                        , "LogIntegralLogFunc   " &
                                                        , "SampleWeight         " &
                                                        , "SampleLogFunc        " &
                                                        ]

#endif

    integer(IK) , parameter :: NUM_DEF_COL = size(COL_HEADER_DEFAULT)   !< the number of columns in the chain file other than the State columns

    type                                    :: ChainFileContents_type
        integer(IK)                         :: ndim = 0_IK
        integer(IK)                         :: lenHeader = 0_IK
        integer(IK)                         :: numDefCol = NUM_DEF_COL
        type(Count_type)                    :: Count
#if defined PARADRAM || defined PARADISE
        integer(IK)         , allocatable   :: Weight(:)                !< The vector of the weights of the MCMC accepted states.
        integer(IK)         , allocatable   :: BurninLoc(:)             !< The burnin locations at the given locations in the chains.
        integer(IK)         , allocatable   :: DelRejStage(:)           !< The delayed rejection stages at which the proposed states were accepted.
        real(RK)            , allocatable   :: Adaptation(:)            !< The vector of the adaptation measures at the MCMC accepted states.
#elif defined PARANEST
        real(RK)            , allocatable   :: LogIntegralLogFunc(:)    !< The natural logarithm of the integral of the user-input target function up to the remaining prior mass.
        real(RK)            , allocatable   :: RemainingPriorMass(:)    !< The remaining prior mass at any stage during the stochastic integration.
        real(RK)            , allocatable   :: Weight(:)                !< The vector of the weights of the MCMC accepted states.
#endif
        real(RK)            , allocatable   :: MeanAccRate(:)           !< The vector of the average acceptance rates at the given point in the chain.
        real(RK)            , allocatable   :: LogFunc(:)               !< The vector of LogFunc values corresponding to the MCMC states.
        real(RK)            , allocatable   :: State(:,:)               !< The (nd,chainSize) MCMC chain of accepted proposed states.
        integer(IK)         , allocatable   :: ProcessID(:)             !< The vector of the ID of the images whose function calls haven been accepted.
        type(CharVec_type)  , allocatable   :: ColHeader(:)             !< The column headers of the chain file.
        character(:)        , allocatable   :: delimiter                !< The delimiter used to separate objects in the chain file.
        type(Err_type)                      :: Err
    contains
        procedure, pass :: nullify => nullifyChainFileContents
        procedure, pass :: get => getChainFileContents
        procedure, pass :: writeChainFile
        procedure, pass :: getLenHeader
        procedure, pass :: writeHeader
    end type ChainFileContents_type

    interface ChainFileContents_type
        module procedure :: constructChainFileContents
    end interface ChainFileContents_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is the constructor of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Return an object of class [ChainFileContents_type](@ref chainfilecontents_type) given the input specifications.
    !>
    !> @param[in]   ndim                : The number of dimensions of the domain of the objective function.
    !> @param[in]   variableNameList    : The list of variable names corresponding to each axis of the domain of the objective function (**optional**).
    !> @param[in]   chainFilePath       : The list of variable names corresponding to each axis of the domain of the objective function (**optional**).
    !> @param[in]   chainSize           : The size of the chain in the chain file specified by the input `chainFilePath` (**optional**).
    !> @param[in]   chainFileForm       : The file format of the chain file (`"binary"` vs. `"compact"` vs. `"verbose"`) (**optional**).
    !> @param[in]   lenHeader           : The full length of the first line in the input file (the header line) (**optional**).
    !> @param[in]   delimiter           : The delimiter symbol used in the chain file (**optional**).
    !> @param[in]   targetChainSize     : The final target size of the chain (in case the chain file is an interrupted simulation) (**optional**).
    !>
    !> \return
    !> `CFC` : An object of class [ChainFileContents_type](@ref chainfilecontents_type) containing the chain.
    !>
    !> \warning
    !> If `chainFilePath` is given, then the rest of the optional arguments *must be also given*.
    function constructChainFileContents(ndim,variableNameList,chainFilePath,chainSize,chainFileForm,lenHeader,delimiter,targetChainSize) result(CFC)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructChainFileContents
#endif
        implicit none
        integer(IK) , intent(in)            :: ndim
        character(*), intent(in), optional  :: chainFileForm
        character(*), intent(in), optional  :: variableNameList(ndim)
        character(*), intent(in), optional  :: chainFilePath
        character(*), intent(in), optional  :: delimiter
        integer(IK) , intent(in), optional  :: lenHeader, chainSize, targetChainSize
        type(ChainFileContents_type)        :: CFC
        type(Err_type)                      :: Err
        integer(IK)                         :: icol
        Err%occurred = .false.

        CFC%ndim = ndim

        ! set up the chain file column header

        allocate(CFC%ColHeader(ndim+NUM_DEF_COL))
        do icol = 1, NUM_DEF_COL
            CFC%ColHeader(icol)%record = trim(adjustl(COL_HEADER_DEFAULT(icol)))
        end do
        if (present(variableNameList)) then
            do icol = NUM_DEF_COL + 1, NUM_DEF_COL + ndim
                CFC%ColHeader(icol)%record = trim(adjustl(variableNameList(icol-NUM_DEF_COL)))
            end do
        end if

        ! set up other variables if given

        if (present(lenHeader)) CFC%lenHeader = lenHeader
        if (present(delimiter)) CFC%delimiter = delimiter
        if (present(targetChainSize)) CFC%Count%target = targetChainSize

        ! read the chain file if the path is given

        if (present(chainFilePath) .and. present(chainFileForm)) call CFC%get(chainFilePath,chainFileForm,Err,chainSize,lenHeader,ndim,delimiter,targetChainSize)
        if (Err%occurred) then
        ! LCOV_EXCL_START
            CFC%Err%occurred = .true.
            CFC%Err%msg = Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

    end function constructChainFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is a method of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Return the contents of a ParaMonte simulation output chain file, always in `compact` format, regardless of the
    !> value of `chainFileFormat` and store it in the object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !>
    !> @param[inout]    CFC             : The object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !> @param[in]       chainFilePath   : The list of variable names corresponding to each axis of the domain of the objective function.
    !> @param[in]       chainFileForm   : The file format of the chain file (`"binary"` vs. `"compact"` vs. `"verbose"`).
    !> @param[out]      Err             : An object of class [Err_type](@ref err_mod::err_type) containing information about whether an error has occurred.
    !> @param[in]       chainSize       : The size of the chain in the chain file specified by the input `chainFilePath` (**optional**).
    !> @param[in]       lenHeader       : The full length of the first line in the input file (the header line) (**optional**).
    !> @param[in]       ndim            : The number of dimensions of the domain of the objective function (**optional**).
    !> @param[in]       delimiter       : The delimiter symbol used in the chain file (**optional**).
    !> @param[in]       targetChainSize : The final target size of the chain (in case the chain file is an interrupted simulation) (**optional**).
    !>
    !> \warning
    !> `targetChainSize` must be `>= chainSize`, if provided. It is used for the allocation of the chain components.
    !> 
    !> \warning
    !> `chainSize` must be `<= targetChainSize`. The first `chainSize` elements of the `CFC` components will contain
    !> the chain information read from the chain file. The chain component elements beyond `chainSize` will be set to zero.
    subroutine getChainFileContents(CFC,chainFilePath,chainFileForm,Err,chainSize,lenHeader,ndim,delimiter,targetChainSize)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getChainFileContents
#endif
        use FileContents_mod, only: getNumRecordInFile
        use Constants_mod, only: IK, RK, NLC, NEGINF_IK, NEGINF_RK
        use String_mod, only: String_type, getLowerCase, num2str

        implicit none

        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME // "@getChainFileContents()"

        class(ChainFileContents_type), intent(inout)    :: CFC
        character(*)    , intent(in)                    :: chainFilePath
        character(*)    , intent(in)                    :: chainFileForm
        type(Err_type)  , intent(out)                   :: Err
        character(*)    , intent(in), optional          :: delimiter
        integer(IK)     , intent(in), optional          :: chainSize, lenHeader, ndim, targetChainSize
        character(:)    , allocatable                   :: chainFilePathTrimmed, thisForm
        character(:)    , allocatable                   :: chainFileFormLowerCase
        type(String_type)                               :: Record
        integer(IK)                                     :: chainFileUnit, i, iState, delimiterLen, chainSizeDefault
        integer(IK)                                     :: irowLastUniqueSample
        integer(IK)                                     :: numColTot
        logical                                         :: fileExists, fileIsOpen, delimHasBegun, delimHasEnded
        logical                                         :: isBinary
        logical                                         :: isCompact
        logical                                         :: isVerbose

        Err%occurred = .false.
        chainFilePathTrimmed = trim(adjustl(chainFilePath))
        inquire(file=chainFilePathTrimmed,exist=fileExists,opened=fileIsOpen,number=chainFileUnit)

        blockFileExistence: if (fileExists) then

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! set up chain file format
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            isBinary = .false.
            isCompact = .false.
            isVerbose = .false.
            chainFileFormLowerCase = getLowerCase(chainFileForm)
            if (chainFileFormLowerCase=="binary") then
                isBinary = .true.
            elseif (chainFileFormLowerCase=="compact" .or. chainFileFormLowerCase=="ascii") then ! "compact" is valid in ParaMCMC, "ascii" is valid in ParaNest
                isCompact = .true.
#if defined PARADRAM || defined PARADISE
            elseif (chainFileFormLowerCase=="verbose") then ! only valid in ParaMCMC
                isVerbose = .true.
#endif
            else
                ! LCOV_EXCL_START
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME//": Unrecognized chain file form: "//chainFileForm
                return
                ! LCOV_EXCL_STOP
            end if

            if (isBinary) then
                thisForm = "unformatted"
                if (.not. present(ndim) .or. .not. present(lenHeader) .or. .not. present(delimiter)) then
                    ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": If the chain file is in binary form, chainSize, lenHeader, delimiter, and ndim must be provided by the user."
                    return
                    ! LCOV_EXCL_STOP
                end if
            else
                thisForm = "formatted"
            end if

            if (fileIsOpen) then
                if (chainFileUnit==-1) then
                    ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": The file located at: "//chainFilePathTrimmed//NLC//"is open, but no unit is connected to the file."//NLC
                    return
                    ! LCOV_EXCL_STOP
                else
                    close(chainFileUnit)
                end if
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! get the number of records in file, minus header line
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (present(chainSize)) then
                chainSizeDefault = chainSize
            else ! here chainSizeDefault is indeed max(chainSize) depending on the file format: verbose or compact
                if (isBinary) then
                    open( newunit = chainFileUnit &
                        , file = chainFilePathTrimmed &
                        , status = "old" &
                        , form = thisForm &
                        , iostat = Err%stat &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                        , SHARED &
#endif
                        )
                    if (Err%stat/=0) then
                        ! LCOV_EXCL_START
                        Err%occurred = .true.
                        Err%msg = PROCEDURE_NAME//": Unable to open the file located at: "//chainFilePathTrimmed//NLC
                        return
                        ! LCOV_EXCL_STOP
                    end if
                    if (allocated(Record%value)) deallocate(Record%value)
                    allocate( character(lenHeader) :: Record%value )
                    read(chainFileUnit) Record%value
                    block
                        integer(IK)             :: processID ! LCOV_EXCL_LINE
#if defined PARADRAM || defined PARADISE
                        integer(IK)             :: delRejStage ! LCOV_EXCL_LINE
                        integer(IK)             :: burninLoc ! LCOV_EXCL_LINE
                        integer(IK)             :: weight ! LCOV_EXCL_LINE
                        real(RK)                :: adaptation ! LCOV_EXCL_LINE
#elif defined PARANEST
                        real(RK)                :: logIntegralLogFunc ! LCOV_EXCL_LINE
                        real(RK)                :: remainingPriorMass ! LCOV_EXCL_LINE
                        real(RK)                :: weight ! LCOV_EXCL_LINE
#endif
                        real(RK)                :: meanAccRate ! LCOV_EXCL_LINE
                        real(RK)                :: logFunc ! LCOV_EXCL_LINE
                        real(RK), allocatable   :: State(:) ! LCOV_EXCL_LINE
                        if (allocated(State)) deallocate(State); allocate(State(ndim))
                        chainSizeDefault = 0_IK
                        loopFindChainSizeDefault: do
#if defined PARADRAM || defined PARADISE
                            read(chainFileUnit,iostat=Err%stat) processID, delRejStage, meanAccRate, adaptation, burninLoc, weight, logFunc, State
#elif defined PARANEST
                            read(chainFileUnit,iostat=Err%stat) processID, meanAccRate, remainingPriorMass, logIntegralLogFunc, weight, logFunc, State
#endif
                            if (Err%stat==0_IK) then
                                chainSizeDefault = chainSizeDefault + 1_IK
                            elseif (is_iostat_end(Err%stat)) then
                                exit loopFindChainSizeDefault
                            ! LCOV_EXCL_START
                            elseif (is_iostat_eor(Err%stat)) then
                                Err%occurred = .true.
                                Err%msg = PROCEDURE_NAME//": Incomplete record detected while reading the input binary chain file at: "//chainFilePathTrimmed//NLC
                                return
                            else
                                Err%occurred = .true.
                                Err%msg = PROCEDURE_NAME//": IO error occurred while reading the input binary chain file at: "//chainFilePathTrimmed//NLC
                                return
                            ! LCOV_EXCL_STOP
                            end if
                        end do loopFindChainSizeDefault
                    end block
                    close(chainFileUnit)
                else ! is not binary
                    call getNumRecordInFile(chainFilePathTrimmed,chainSizeDefault,Err,exclude="")
                    if (Err%occurred) then
                    ! LCOV_EXCL_START
                        Err%msg = PROCEDURE_NAME//Err%msg
                        return
                    end if
                    ! LCOV_EXCL_STOP
                    chainSizeDefault = chainSizeDefault - 1_IK ! subtract header
                end if
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! set the number of elements in the Chain components
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (present(targetChainSize)) then ! in restart mode, this must always be the case
                CFC%Count%target = targetChainSize
            else
                CFC%Count%target = chainSizeDefault
            end if
            !if (CFC%Count%target<chainSizeDefault) then
            !    Err%occurred = .true.
            !    Err%msg =   PROCEDURE_NAME//": Internal error occurred. The input targetChainSize cannot be smaller than the input chainSize:" // NLC // &
            !                "    targetChainSize = " // num2str(CFC%Count%target) // NLC // &
            !                "          chainSize = " // num2str(chainSizeDefault) // NLC // &
            !                "It appears that the user has manipulated the output chain file."
            !    return
            !end if

            ! allocate Chain components

#if defined PARADRAM || defined PARADISE
            if (allocated(CFC%DelRejStage))         deallocate(CFC%DelRejStage);        allocate(CFC%DelRejStage        (CFC%Count%target), source = NEGINF_IK)
            if (allocated(CFC%Adaptation))          deallocate(CFC%Adaptation);         allocate(CFC%Adaptation         (CFC%Count%target), source = NEGINF_RK) ! this initialization is critical and relied upon later below
            if (allocated(CFC%BurninLoc))           deallocate(CFC%BurninLoc);          allocate(CFC%BurninLoc          (CFC%Count%target), source = NEGINF_IK)
#elif defined PARANEST
            if (allocated(CFC%RemainingPriorMass))  deallocate(CFC%RemainingPriorMass)  allocate(CFC%RemainingPriorMass (CFC%Count%target), source = NEGINF_RK)
            if (allocated(CFC%LogIntegralLogFunc))  deallocate(CFC%LogIntegralLogFunc)  allocate(CFC%LogIntegralLogFunc (CFC%Count%target), source = NEGINF_RK)
#endif
            if (allocated(CFC%MeanAccRate))         deallocate(CFC%MeanAccRate);        allocate(CFC%MeanAccRate        (CFC%Count%target), source = NEGINF_RK)
            if (allocated(CFC%ProcessID))           deallocate(CFC%ProcessID);          allocate(CFC%ProcessID          (CFC%Count%target), source = NEGINF_IK)
            if (allocated(CFC%LogFunc))             deallocate(CFC%LogFunc);            allocate(CFC%LogFunc            (CFC%Count%target), source = NEGINF_RK)
            if (allocated(CFC%Weight))              deallocate(CFC%Weight);             allocate(CFC%Weight             (CFC%Count%target), source = NEGINF_IK)
            if (allocated(CFC%State))               deallocate(CFC%State)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! find the delimiter
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockFindDelim: if (present(delimiter)) then

                CFC%delimiter = delimiter

            else blockFindDelim

                if (allocated(CFC%delimiter)) deallocate(CFC%delimiter)
                allocate( character(1023) :: CFC%delimiter )
                if (allocated(Record%value)) deallocate(Record%value)
                allocate( character(99999) :: Record%value )

                open( newunit = chainFileUnit &
                    , file = chainFilePathTrimmed &
                    , status = "old" &
                    , form = thisForm &
                    , iostat = Err%stat &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                    , SHARED &
#endif
                    )
                if (Err%stat/=0) then
                ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": Unable to open the file located at: "//chainFilePathTrimmed//"."//NLC
                    return
                end if
                ! LCOV_EXCL_STOP

                read(chainFileUnit,*)   ! skip the header
                read(chainFileUnit,"(A)") Record%value  ! read the first numeric row in string format
                close(chainFileUnit)

                Record%value = trim(adjustl(Record%value))
                delimHasEnded = .false.
                delimHasBegun = .false.
                delimiterLen = 0
                loopSearchDelimiter: do i = 1, len(Record%value)-1
                    if ( Record%isDigit(Record%value(i:i)) ) then
                        if (delimHasBegun) delimHasEnded = .true.
                    elseif (Record%value(i:i)=="." .or. Record%value(i:i)=="+" .or. Record%value(i:i)=="-") then
                        if (delimHasBegun) then
                            delimHasEnded = .true.
                        else
                            ! LCOV_EXCL_START
                            Err%occurred = .true.
                            Err%msg = PROCEDURE_NAME//": The file located at: " // chainFilePathTrimmed //NLC//&
                            "has unrecognizable format. Found "//Record%value(i:i)//" in the first column, while expecting positive integer."//NLC
                            return
                            ! LCOV_EXCL_STOP
                        end if
                    else
                        if (i==1) then  ! here it is assumed that the first column in chain file always contains integers
                            ! LCOV_EXCL_START
                            Err%occurred = .true.
                            Err%msg = PROCEDURE_NAME//": The file located at: "//chainFilePathTrimmed//NLC//"has unrecognizable format."//NLC
                            return
                            ! LCOV_EXCL_STOP
                        else
                            delimHasBegun = .true.
                            delimiterLen = delimiterLen + 1
                            CFC%delimiter(delimiterLen:delimiterLen) = Record%value(i:i)
                        end if
                    end if
                    if (delimHasEnded) exit loopSearchDelimiter
                end do loopSearchDelimiter

                if (.not.(delimHasBegun.and.delimHasEnded)) then
                    ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": The file located at: "//chainFilePathTrimmed//NLC//"has unrecognizable format. Could not identify the column delimiter."//NLC
                    return
                    ! LCOV_EXCL_STOP
                else
                    CFC%delimiter = trim(adjustl(CFC%delimiter(1:delimiterLen)))
                    delimiterLen = len(CFC%delimiter)
                    if (delimiterLen==0) then
                        CFC%delimiter = " "
                        delimiterLen = 1
                    end if
                end if

            end if blockFindDelim

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! find the number of dimensions of the state (the number of function variables)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (present(ndim)) then
                CFC%ndim = ndim
            else
                Record%Parts = Record%split(Record%value,CFC%delimiter,Record%nPart)
                CFC%numDefCol = 0_IK
                loopFindNumDefCol: do i = 1, Record%nPart
                    if ( index(string=Record%Parts(i)%record,substring="LogFunc") > 0 ) then
                        CFC%numDefCol = i
                        exit loopFindNumDefCol
                    end if
                end do loopFindNumDefCol
                if (CFC%numDefCol/=NUM_DEF_COL .or. CFC%numDefCol==0_IK) then
                    ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": Internal error occurred. CFC%numDefCol/=NUM_DEF_COL: " // num2str(CFC%numDefCol) // num2str(NUM_DEF_COL)
                    return
                    ! LCOV_EXCL_STOP
                end if
                CFC%ndim = Record%nPart - NUM_DEF_COL
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! reopen the file to read the contents
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            open( newunit = chainFileUnit &
                , file = chainFilePathTrimmed &
                , status = "old" &
                , form = thisForm &
                , iostat = Err%stat &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                , SHARED &
#endif
                )
            if (Err%stat/=0) then
                ! LCOV_EXCL_START
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME//": Unable to open the file located at: "//chainFilePathTrimmed //"."//NLC
                return
                ! LCOV_EXCL_STOP
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! first read the column headers
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (allocated(Record%value)) deallocate(Record%value) ! set up the record string that keeps the contents of each line
            if (isBinary) then
                allocate( character(lenHeader) :: Record%value )
                read(chainFileUnit) Record%value
            else
                allocate( character(99999) :: Record%value ) ! such huge allocation is rather redundant and is good for a ~4000 dimensional objective function.
                read(chainFileUnit, "(A)" ) Record%value
            end if
            CFC%ColHeader = Record%split(trim(adjustl(Record%value)), CFC%delimiter, Record%npart)
            do i = 1, Record%npart ! xxx is this trimming necessary?
                CFC%ColHeader(i)%record = trim(adjustl(CFC%ColHeader(i)%record))
            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! read the chain
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (.not. isBinary) then
                numColTot = CFC%numDefCol + CFC%ndim
            end if

            CFC%Count%verbose = 0._RK
            allocate(CFC%State(CFC%ndim,CFC%Count%target))
            blockChainFileFormat: if (isBinary) then

                loopReadBinary: do iState = 1, chainSizeDefault
#if defined PARADRAM || defined PARADISE
                    read(chainFileUnit, iostat=Err%stat ) CFC%ProcessID                (iState)    &
                                                        , CFC%DelRejStage              (iState)    &
                                                        , CFC%MeanAccRate              (iState)    &
                                                        , CFC%Adaptation               (iState)    &
                                                        , CFC%BurninLoc                (iState)    &
                                                        , CFC%Weight                   (iState)    &
                                                        , CFC%LogFunc                  (iState)    &
                                                        , CFC%State         (1:CFC%ndim,iState)
#elif defined PARANEST
                    read(chainFileUnit, iostat=Err%stat ) CFC%ProcessID                (iState)    &
                                                        , CFC%MeanAccRate              (iState)    &
                                                        , CFC%RemainingPriorMass       (iState)    &
                                                        , CFC%LogIntegralLogFunc       (iState)    &
                                                        , CFC%Weight                   (iState)    &
                                                        , CFC%LogFunc                  (iState)    &
                                                        , CFC%State         (1:CFC%ndim,iState)
#endif
                    if (is_iostat_eor(Err%stat) .or. is_iostat_end(Err%stat)) then
                    ! LCOV_EXCL_START
                        call warnUserAboutCorruptChainFile(iState)
                        exit loopReadBinary
                    end if
                    ! LCOV_EXCL_STOP
                    CFC%Count%verbose = CFC%Count%verbose + CFC%Weight(iState)
                end do loopReadBinary

            elseif (isCompact) then blockChainFileFormat

                loopReadCompact: do iState = 1, chainSizeDefault
                    read(chainFileUnit, "(A)" ) Record%value
                    Record%Parts = Record%split(trim(adjustl(Record%value)),CFC%delimiter,Record%nPart)
                    if (Record%nPart<numColTot) then
                        ! LCOV_EXCL_START
                        call warnUserAboutCorruptChainFile(iState)
                        exit loopReadCompact
                        ! LCOV_EXCL_STOP
                    else
#if defined PARADRAM || defined PARADISE
                        read(Record%Parts(1)%record,*) CFC%ProcessID            (iState)
                        read(Record%Parts(2)%record,*) CFC%DelRejStage          (iState)
                        read(Record%Parts(3)%record,*) CFC%MeanAccRate          (iState)
                        read(Record%Parts(4)%record,*) CFC%Adaptation           (iState)
                        read(Record%Parts(5)%record,*) CFC%BurninLoc            (iState)
                        read(Record%Parts(6)%record,*) CFC%Weight               (iState)
                        read(Record%Parts(7)%record,*) CFC%LogFunc              (iState)
#elif defined PARANEST
                        read(Record%Parts(1)%record,*) CFC%ProcessID            (iState)
                        read(Record%Parts(2)%record,*) CFC%MeanAccRate          (iState)
                        read(Record%Parts(3)%record,*) CFC%RemainingPriorMass   (iState)
                        read(Record%Parts(4)%record,*) CFC%LogIntegralLogFunc   (iState)
                        read(Record%Parts(5)%record,*) CFC%Weight               (iState)
                        read(Record%Parts(6)%record,*) CFC%LogFunc              (iState)
#endif
                        do i = 1, CFC%ndim
                            read(Record%Parts(NUM_DEF_COL+i)%record,*) CFC%State(i,iState)
                        end do
                        CFC%Count%verbose = CFC%Count%verbose + CFC%Weight(iState)
                    end if
                end do loopReadCompact

#if defined PARADRAM || defined PARADISE

            else blockChainFileFormat ! is verbose form

                blockChainSizeDefault: if (chainSizeDefault>0_IK) then

                    CFC%Count%compact = 1_IK
                    blockReadVerbose: block

                        logical                 :: newUniqueSampleDetected
                        integer(IK)             :: processID
                        integer(IK)             :: delRejStage
                        real(RK)                :: meanAccRate
                        real(RK)                :: adaptation
                        integer(IK)             :: burninLoc
                        integer(IK)             :: weight
                        real(RK)                :: logFunc
                        real(RK), allocatable   :: State(:)
                        if (allocated(State)) deallocate(State); allocate(State(ndim))

                        irowLastUniqueSample = 0_IK

                        ! read the first sample

                        read(chainFileUnit, "(A)" ) Record%value
                        Record%Parts = Record%split(trim(adjustl(Record%value)),CFC%delimiter,Record%nPart)
                        if (Record%nPart<numColTot) then
                            ! LCOV_EXCL_START
                            call warnUserAboutCorruptChainFile(iState)
                            !exit blockChainSizeDefault
                            ! intel 2018 to 2019.05 yields internal compiler error with the above exit. Intel 19.1 and gnu 9.1 are fine.
                            ! The following is a workaround for now.
                            exit blockReadVerbose
                            ! LCOV_EXCL_STOP
                        else
                            read(Record%Parts(1)%record,*) CFC%ProcessID(CFC%Count%compact)
                            read(Record%Parts(2)%record,*) CFC%DelRejStage(CFC%Count%compact)
                            read(Record%Parts(3)%record,*) CFC%MeanAccRate(CFC%Count%compact)
                            read(Record%Parts(4)%record,*) CFC%Adaptation(CFC%Count%compact)
                            read(Record%Parts(5)%record,*) CFC%BurninLoc(CFC%Count%compact)
                            read(Record%Parts(6)%record,*) CFC%Weight(CFC%Count%compact)
                            read(Record%Parts(7)%record,*) CFC%LogFunc(CFC%Count%compact)
                            do i = 1, CFC%ndim
                                read(Record%Parts(CFC%numDefCol+i)%record,*) CFC%State(i,CFC%Count%compact)
                            end do
                        end if

                        ! read the rest of samples beyond the first, if any exist

                        newUniqueSampleDetected = .false.
                        loopOverChainfFileContents: do iState = 2, chainSizeDefault

                            read(chainFileUnit, "(A)" ) Record%value
                            Record%Parts = Record%split(trim(adjustl(Record%value)), CFC%delimiter, Record%nPart)
                            if (Record%nPart<numColTot) then
                                ! LCOV_EXCL_START
                                call warnUserAboutCorruptChainFile(iState)
                                exit loopOverChainfFileContents
                                ! LCOV_EXCL_STOP
                            else
                                read(Record%Parts(1)%record,*) ProcessID
                                read(Record%Parts(2)%record,*) DelRejStage
                                read(Record%Parts(3)%record,*) MeanAccRate
                                read(Record%Parts(4)%record,*) Adaptation
                                read(Record%Parts(5)%record,*) BurninLoc
                                read(Record%Parts(6)%record,*) Weight
                                read(Record%Parts(7)%record,*) LogFunc
                                do i = 1, CFC%ndim
                                    read(Record%Parts(CFC%numDefCol+i)%record,*) State(i)
                                end do

                                ! increment CFC%Count%compact if new sample detected

                                newUniqueSampleDetected =    LogFunc        /= CFC%LogFunc    (CFC%Count%compact) &
                                                       !.or. MeanAccRate    /= CFC%MeanAccRate(CFC%Count%compact) &
                                                       !.or. Adaptation     /= CFC%Adaptation (CFC%Count%compact) &
                                                       !.or. BurninLoc      /= CFC%BurninLoc  (CFC%Count%compact) &
                                                       !.or. Weight         /= CFC%Weight     (CFC%Count%compact) &
                                                       !.or. DelRejStage    /= CFC%DelRejStage(CFC%Count%compact) &
                                                       !.or. ProcessID      /= CFC%ProcessID  (CFC%Count%compact) &
                                                        .or. any(CFC%State(1:CFC%ndim,CFC%Count%compact) /= CFC%State(1:CFC%ndim,CFC%Count%compact))
                                if (newUniqueSampleDetected) then
                                    irowLastUniqueSample = irowLastUniqueSample + CFC%Weight(CFC%Count%compact)
                                    ! increment the compact sample
                                    CFC%Count%compact = CFC%Count%compact + 1_IK
                                    if (CFC%Count%target<CFC%Count%compact) then
                                        Err%occurred = .true.
                                        Err%msg =   PROCEDURE_NAME//": Fatal error occurred. CFC%Count%target<CFC%Count%compact: "// &
                                                    num2str(CFC%Count%target) // " /= " // num2str(CFC%Count%compact) // &
                                                    "The contents of the input chain file is longer than the user-requested allocation size."
                                        return
                                    end if
                                else
                                    weight = CFC%Weight(CFC%Count%compact) + 1_IK
                                end if

                                ! write the latest sample

                                CFC%LogFunc         (CFC%Count%compact) = LogFunc
                                CFC%MeanAccRate     (CFC%Count%compact) = MeanAccRate
                                CFC%Adaptation      (CFC%Count%compact) = max(CFC%Adaptation(CFC%Count%compact),Adaptation)
                                CFC%BurninLoc       (CFC%Count%compact) = BurninLoc
                                CFC%Weight          (CFC%Count%compact) = Weight
                                CFC%DelRejStage     (CFC%Count%compact) = DelRejStage
                                CFC%ProcessID       (CFC%Count%compact) = ProcessID
                                CFC%State(1:CFC%ndim,CFC%Count%compact) = State(1:CFC%ndim)

                            end if

                        end do loopOverChainfFileContents

                    end block blockReadVerbose

                else blockChainSizeDefault

                    CFC%Count%compact = 0_IK
                    CFC%Count%verbose = 0_IK

                end if blockChainSizeDefault
#endif

            end if blockChainFileFormat

            if (isBinary .or. isCompact) then
                CFC%Count%compact = chainSizeDefault
            else
                CFC%Count%verbose = chainSizeDefault
                if (CFC%Count%verbose/=sum(CFC%Weight(1:CFC%Count%compact))) then
                    ! LCOV_EXCL_START
                    Err%occurred = .true.
                    Err%msg =   PROCEDURE_NAME//": Internal error occurred. CountVerbose/=sum(Weight): "// &
                                num2str(CFC%Count%verbose)//" /= "//num2str(sum(CFC%Weight(1:CFC%Count%compact)))// &
                                ", CFC%Count%compact = "//num2str(CFC%Count%compact)
                    return
                    ! LCOV_EXCL_STOP
                elseif (.not. present(targetChainSize)) then
                    CFC%ProcessID     = CFC%ProcessID   (1:CFC%Count%compact)
                    CFC%DelRejStage   = CFC%DelRejStage (1:CFC%Count%compact)
                    CFC%MeanAccRate   = CFC%MeanAccRate (1:CFC%Count%compact)
                    CFC%Adaptation    = CFC%Adaptation  (1:CFC%Count%compact)
                    CFC%BurninLoc     = CFC%BurninLoc   (1:CFC%Count%compact)
                    CFC%Weight        = CFC%Weight      (1:CFC%Count%compact)
                    CFC%LogFunc       = CFC%LogFunc     (1:CFC%Count%compact)
                    CFC%State         = CFC%State       (1:CFC%ndim,1:CFC%Count%compact)
                end if
            end if

            close(chainFileUnit)

            ! set the rest of elements to null values

            if (CFC%Count%target>chainSizeDefault) call CFC%nullify(startIndex=CFC%Count%compact+1_IK, endIndex=CFC%Count%target)

        else blockFileExistence

            ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME//": The chain file does not exist in the given file path: "//chainFilePathTrimmed
            return
            ! LCOV_EXCL_STOP

        end if blockFileExistence

    contains

        ! LCOV_EXCL_START
        subroutine warnUserAboutCorruptChainFile(lineNumber)
            implicit none
            integer(IK) :: lineNumber
            if (isVerbose) then
                chainSizeDefault = irowLastUniqueSample
                CFC%Count%compact = CFC%Count%compact - 1
            else
                chainSizeDefault = chainSizeDefault - 1
            end if
            call warn   ( prefix = INDENT//"ParaMonte" &
                        , marginTop = 0_IK &
                        , marginBot = 2_IK &
                        , outputUnit = output_unit &
                        , msg = "An end-of-file or end-of-record condition occurred while parsing the contents of the chain file at line = "//num2str(lineNumber)//" with iostat = "//num2str(Err%stat)// &
                                ". Assuming the previous line as the last line of the chain file..." &
                        )
        ! LCOV_EXCL_STOP
        end subroutine warnUserAboutCorruptChainFile

    end subroutine getChainFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is a method of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Reset the components of the chain object to an unlikely value for the purpose of error catching and debugging.
    !> Store the modified components as part of the input object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !>
    !> @param[inout]    CFC             : The number of dimensions of the domain of the objective function.
    !> @param[in]       startIndex      : The beginning index beyond which the component values will be reset.
    !> @param[in]       endIndex        : The ending index below which the component values will be reset.
    subroutine nullifyChainFileContents(CFC,startIndex,endIndex)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyChainFileContents
#endif
        implicit none
        class(ChainFileContents_type), intent(inout)    :: CFC
        integer(IK), intent(in)                         :: startIndex, endIndex
#if defined PARADRAM || defined PARADISE
        CFC%ProcessID           (startIndex:endIndex) = -huge(0_IK)
        CFC%DelRejStage         (startIndex:endIndex) = -huge(0_IK)
        CFC%MeanAccRate         (startIndex:endIndex) = -huge(0._RK)
        CFC%Adaptation          (startIndex:endIndex) = -huge(0._RK)
        CFC%BurninLoc           (startIndex:endIndex) = -huge(0_IK)
        CFC%Weight              (startIndex:endIndex) = 0_IK
        CFC%LogFunc             (startIndex:endIndex) = -huge(0._RK)
        CFC%State               (1:CFC%ndim,startIndex:endIndex) = -huge(0._RK)
#elif defined PARANEST
        CFC%ProcessID           (startIndex:endIndex) = -huge(0_IK)
        CFC%MeanAccRate         (startIndex:endIndex) = -huge(0._RK)
        CFC%RemainingPriorMass  (startIndex:endIndex) = -huge(0._RK)
        CFC%LogIntegralLogFunc  (startIndex:endIndex) = -huge(0_IK)
        CFC%Weight              (startIndex:endIndex) = 0._RK
        CFC%LogFunc             (startIndex:endIndex) = -huge(0._RK)
        CFC%State               (1:CFC%ndim,startIndex:endIndex) = -huge(0._RK)
#endif
    end subroutine nullifyChainFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is a method of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Return the length of the header of the chain file.
    !>
    !> @param[inout]    CFC             :   The object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !> @param[in]       ndim            :   The number of dimensions of the domain of the objective function.
    !> @param[in]       isBinary        :   The logical flag indicating whether the file is in `binary` format.
    !> @param[in]       chainFileFormat :   The Fortran IO formatting string to be used to read the contents of the chain file (**optional**).
    !>                                      This argument is only required with a non-binary chain file, i.e., when `isBinary = .false.`.
    subroutine getLenHeader(CFC,ndim,isBinary,chainFileFormat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenHeader
#endif
        use Constants_mod, only: IK ! LCOV_EXCL_LINE
        use Err_mod, only: abort
        implicit none
        class(ChainFileContents_type), intent(inout)    :: CFC
        integer(IK) , intent(in)                        :: ndim
        logical     , intent(in)                        :: isBinary
        character(*), intent(in), optional              :: chainFileFormat
        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME//"@getLenHeader()"
        character(:), allocatable                       :: record
        integer(IK)                                     :: i
        CFC%Err%occurred = .false.
        allocate( character(99999) :: record )
        if (isBinary) then
            write( record , "(*(g0,:,','))" ) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
        else
            if ( present(chainFileFormat) ) then
                write(record,chainFileFormat) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
            else
                ! LCOV_EXCL_START
                CFC%Err%occurred = .true.
                CFC%Err%msg = PROCEDURE_NAME//"Internal error occurred. For formatted chain files, chainFileFormat must be given."
                call abort(CFC%Err)
                error stop
                return
                ! LCOV_EXCL_STOP
            end if
        end if
        CFC%lenHeader = len_trim(adjustl(record))
        deallocate(record)
    end subroutine getLenHeader

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is a method of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Write the requested header to the chain file.
    !>
    !> @param[inout]    CFC             :   The object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !> @param[in]       ndim            :   The number of dimensions of the domain of the objective function.
    !> @param[in]       chainFileUnit   :   The unit ID of the chain file to which the header should be written.
    !> @param[in]       isBinary        :   The logical flag indicating whether the file is in `binary` format.
    !> @param[in]       chainFileFormat :   The Fortran IO formatting string to be used to read the contents of the chain file (**optional**).
    !>                                      This argument is only required with a non-binary chain file, i.e., when `isBinary = .false.`.
    subroutine writeHeader(CFC,ndim,chainFileUnit,isBinary,chainFileFormat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writeHeader
#endif
        use Constants_mod, only: IK
        use Err_mod, only: abort
        implicit none
        class(ChainFileContents_type), intent(inout)    :: CFC
        integer(IK) , intent(in)                        :: ndim, chainFileUnit
        logical     , intent(in)                        :: isBinary
        character(*), intent(in), optional              :: chainFileFormat
        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME//"@writeHeader()"
        character(:), allocatable                       :: record
        integer(IK)                                     :: i
        CFC%Err%occurred = .false.
        if (isBinary) then
            allocate( character(99999) :: record )
            write( record , "(*(g0,:,','))" ) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
            write(chainFileUnit) trim(adjustl(record))
            deallocate(record)
        else
            if ( present(chainFileFormat) ) then
                write(chainFileUnit,chainFileFormat) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
            else
                ! LCOV_EXCL_START
                CFC%Err%occurred = .true.
                CFC%Err%msg = PROCEDURE_NAME//"Internal error occurred. For formatted chain files, chainFileFormat must be given."
                call abort(CFC%Err)
                error stop
                return
                ! LCOV_EXCL_STOP
            end if
        end if
    end subroutine writeHeader

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is a method of the class [ChainFileContents_type](@ref chainfilecontents_type).\n
    !> Write the chain properties to the chain file.
    !>
    !> @param[inout]    CFC                     :   The object of class [ChainFileContents_type](@ref chainfilecontents_type).
    !> @param[in]       ndim                    :   The number of dimensions of the domain of the objective function.
    !> @param[in]       compactStartIndex       :   The beginning index of the compact chain beyond which the elements of the chain will be written to the output file.
    !> @param[in]       compactEndIndex         :   The ending index of the compact chain below which the elements of the chain will be written to the output file.
    !> @param[in]       chainFileUnit           :   The unit ID of the chain file to which the header should be written.
    !> @param[in]       chainFileForm           :   The file format of the chain file (`"binary"` vs. `"compact"` vs. `"verbose"`).
    !> @param[in]       chainFileFormat         :   The Fortran IO formatting string to be used to read the contents of the chain file (**optional**).
    !>                                              This argument is only required with a non-binary chain file, i.e., when `isBinary = .false.`.
    !> @param[in]       adaptiveUpdatePeriod    :   The adaptive update period (**optional**). It must be provided if `chainFileForm = "verbose"`.
    subroutine writeChainFile(CFC,ndim,compactStartIndex,compactEndIndex,chainFileUnit,chainFileForm,chainFileFormat,adaptiveUpdatePeriod)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writeChainFile
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: abort
        implicit none
        class(ChainFileContents_type), intent(inout)    :: CFC
        integer(IK) , intent(in)                        :: ndim, compactStartIndex, compactEndIndex, chainFileUnit
        character(*), intent(in)                        :: chainFileForm
        character(*), intent(in), optional              :: chainFileFormat
        integer(IK) , intent(in), optional              :: adaptiveUpdatePeriod
        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME//"@writeChainFile()"
        logical                                         :: isBinary, isCompact, isVerbose
        real(RK)                                        :: adaptation
        integer(IK)                                     :: i,j, counter

        CFC%Err%occurred = .false.

        isBinary = .false.
        isCompact = .false.
        isVerbose = .false.
        if (chainFileForm=="binary") then
            isBinary = .true.
        else
            if (present(chainFileFormat)) then
                if (chainFileForm=="compact") then
                    isCompact = .true.
#if defined PARADRAM || defined PARADISE
                elseif (chainFileForm=="verbose") then
                    if (present(adaptiveUpdatePeriod)) then
                        isVerbose = .true.
                    else
                        ! LCOV_EXCL_START
                        CFC%Err%occurred = .true.
                        CFC%Err%msg = PROCEDURE_NAME//"Internal error occurred. For verbose chain files, adaptiveUpdatePeriod must be given."
                        ! LCOV_EXCL_STOP
                    end if
#endif
                else
                    ! LCOV_EXCL_START
                    CFC%Err%occurred = .true.
                    CFC%Err%msg = PROCEDURE_NAME//"Internal error occurred. Unknown chain file format: "//chainFileForm
                    ! LCOV_EXCL_STOP
                end if
            else
                ! LCOV_EXCL_START
                CFC%Err%occurred = .true.
                CFC%Err%msg = PROCEDURE_NAME//"Internal error occurred. For formatted chain files, chainFileFormat must be given."
                ! LCOV_EXCL_STOP
            end if
        end if

        if (CFC%Err%occurred) then
            ! LCOV_EXCL_START
            call abort(CFC%Err)
            return
            ! LCOV_EXCL_STOP
        end if

        call CFC%writeHeader(ndim,chainFileUnit,isBinary,chainFileFormat)

        if (compactStartIndex<=compactEndIndex) then
            blockChainFileFormat: if (isCompact) then
#if defined PARANEST
                do i = compactStartIndex, compactEndIndex
                    write(chainFileUnit,chainFileFormat     ) CFC%ProcessID(i)          &
                                                            , CFC%MeanAccRate(i)        &
                                                            , CFC%RemainingPriorMass(i) &
                                                            , CFC%LogIntegralLogFunc(i) &
                                                            , CFC%Weight(i)             &
                                                            , CFC%LogFunc(i)            &
                                                            , CFC%State(1:ndim,i)
                end do
#elif defined PARADRAM || defined PARADISE
                do i = compactStartIndex, compactEndIndex
                    write(chainFileUnit,chainFileFormat     ) CFC%ProcessID(i)      &
                                                            , CFC%DelRejStage(i)    &
                                                            , CFC%MeanAccRate(i)    &
                                                            , CFC%Adaptation(i)     &
                                                            , CFC%BurninLoc(i)      &
                                                            , CFC%Weight(i)         &
                                                            , CFC%LogFunc(i)        &
                                                            , CFC%State(1:ndim,i)
                end do
#endif
            elseif (isBinary) then blockChainFileFormat
#if defined PARANEST
                do i = compactStartIndex, compactEndIndex
                    write(chainFileUnit                     ) CFC%ProcessID(i)          &
                                                            , CFC%MeanAccRate(i)        &
                                                            , CFC%RemainingPriorMass(i) &
                                                            , CFC%LogIntegralLogFunc(i) &
                                                            , CFC%Weight(i)             &
                                                            , CFC%LogFunc(i)            &
                                                            , CFC%State(1:ndim,i)
                end do
#elif defined PARADRAM || defined PARADISE
                do i = compactStartIndex, compactEndIndex
                    write(chainFileUnit                     ) CFC%ProcessID(i)      &
                                                            , CFC%DelRejStage(i)    &
                                                            , CFC%MeanAccRate(i)    &
                                                            , CFC%Adaptation(i)     &
                                                            , CFC%BurninLoc(i)      &
                                                            , CFC%Weight(i)         &
                                                            , CFC%LogFunc(i)        &
                                                            , CFC%State(1:ndim,i)
                end do
            elseif (isVerbose) then blockChainFileFormat
                counter = compactStartIndex
                do i = compactStartIndex, compactEndIndex
                    do j = 1, CFC%Weight(i)
                        if (mod(counter,adaptiveUpdatePeriod)==0_IK) then
                            adaptation = CFC%Adaptation(i)
                        else
                            adaptation = 0._RK
                        end if
                        write(chainFileUnit,chainFileFormat ) CFC%ProcessID(i)      &
                                                            , CFC%DelRejStage(i)    &
                                                            , CFC%MeanAccRate(i)    &
                                                            , adaptation            &
                                                            , CFC%BurninLoc(i)      &
                                                            , 1_IK                  &
                                                            , CFC%LogFunc(i)        &
                                                            , CFC%State(1:ndim,i)
                        counter = counter + 1
                    end do
                end do
#endif
            end if blockChainFileFormat
        end if
        flush(chainFileUnit)
    end subroutine writeChainFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaXXXX
