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

module ParaMCMCRefinedChain_mod

    use SpecMCMC_SampleRefinementMethod_mod, only: BATCH_MEANS_METHOD_NAME
    use SpecMCMC_SampleRefinementMethod_mod, only: CUTOFF_AUTOCORR_METHOD_NAME
    use SpecMCMC_SampleRefinementMethod_mod, only: MAX_CUMSUM_AUTOCORR_METHOD_NAME
    use ParaMonteChainFileContents_mod, only: Count_type
    use JaggedArray_mod, only: CharVec_type
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaMCMCRefinedChain_mod"

    type                                    :: RefinedChain_type
        integer(IK)                         :: ndim  = 0_IK         ! number of sampling variables
        integer(IK)                         :: numRefinement = 0_IK ! number of refinements, zero if sample size is prescribed by the user
        type(Count_type)    , allocatable   :: Count(:)             ! compact and verbose counts
        real(RK)            , allocatable   :: IAC(:,:)             ! size of (ndim,0:numRefinement): The Integrated AutoCorrelation Time at each refinement stage
        real(RK)            , allocatable   :: LogFuncState(:,:)    ! size of (ndim,Count%compact): LogFuncState is LogFunc + Variables
        integer(IK)         , allocatable   :: Weight(:)            ! size of (Count%compact): Weight of each state
        type(CharVec_type)  , allocatable   :: ColHeader(:)         ! refined sample column headers
        type(Err_type)                      :: Err
    contains
        procedure, pass :: get => getRefinedChain
        procedure, pass :: write => writeRefinedChain
    end type RefinedChain_type

   type Method_type
       logical :: isBatchMeans = .false.
       logical :: isCutoffAutoCorr = .false.
       logical :: isViaCompactChain = .true.
       logical :: isViaVerboseChain = .true.
       logical :: isMaxCumSumAutoCorr = .false.
   end type Method_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getRefinedChain  ( RefinedChain              &
                                , CFC                       &
                                , Err                       &
                                , burninLoc                 &
                                , refinedChainSize          &
                                , sampleRefinementCount     &
                                , sampleRefinementMethod    &
                                )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefinedChain
#endif

        use, intrinsic :: iso_fortran_env, only: output_unit
        use ParaMonteChainFileContents_mod, only: ChainFileContents_type, NUM_DEF_COL
        use CrossCorr_mod, only: getBatchMeansIAC, getCumSumIAC, getMaxCumSumIAC
        use String_mod, only: getLowerCase

        implicit none

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@getRefinedChain()"

        class(RefinedChain_type)    , intent(inout)             :: RefinedChain
        type(ChainFileContents_type), intent(inout)             :: CFC
        type(Err_type)              , intent(out)               :: Err
        character(*)                , intent(in), optional      :: sampleRefinementMethod
        integer(IK)                 , intent(in), optional      :: burninLoc, refinedChainSize, sampleRefinementCount
        integer(IK)                                             :: burninLocDefault, i, countCompact
        integer(IK)                                             :: maxRefinementCurrent, maxRefinementCount
        real(RK)                                                :: integratedAutoCorrTime
        logical                                                 :: refinedChainSizeIsPresent
        logical                                                 :: maxRefinementCountIsReached
        character(:)    , allocatable                           :: sampleRefinementMethodLowerCase
        type(Count_type), allocatable                           :: DumCount(:)
        real(RK)        , allocatable                           :: DumIAC(:,:), DumVec(:)
        integer(IK)     , allocatable                           :: SampleWeight(:)
        type(Method_type)                                       :: Method
        type(ChainFileContents_type)                            :: DumCFC

        Err%occurred = .false.

        ! if the size of the refined sample is given as input, then generate the requested sample straight

        refinedChainSizeIsPresent = present(refinedChainSize)
        if (refinedChainSizeIsPresent) then    ! ignore sampleRefinementCount, even if it is given by the user
            maxRefinementCount = 1_IK
        else
            maxRefinementCount = 20_IK  ! this is a temporary maximum value, to be increased later if needed
            if (present(sampleRefinementCount)) maxRefinementCount = sampleRefinementCount
        end if

        ! this is to avoid memory overflow due to extremely large maxRefinementCount requested by the user

        maxRefinementCurrent = min(2_IK, maxRefinementCount)

        ! compute ndim and the initial chain size

        RefinedChain%numRefinement = 0_IK

        if (CFC%ndim == 0_IK) then
            RefinedChain%ndim = size(CFC%State(:,1))
            CFC%ndim = RefinedChain%ndim
        else
            RefinedChain%ndim = CFC%ndim
        end if

        ! allocate components

        if (allocated(RefinedChain%IAC)) deallocate(RefinedChain%IAC)
        if (allocated(RefinedChain%Count)) deallocate(RefinedChain%Count)
        allocate(RefinedChain%IAC(0:RefinedChain%ndim,0:maxRefinementCurrent))
        allocate(RefinedChain%Count(0:maxRefinementCurrent))
        if (CFC%Count%compact == 0_IK) then
            RefinedChain%Count(RefinedChain%numRefinement)%compact = size(CFC%State(1,:))
            CFC%Count%compact = RefinedChain%Count(RefinedChain%numRefinement)%compact
        else
            RefinedChain%Count(RefinedChain%numRefinement)%compact = CFC%Count%compact
        end if
        if (CFC%Count%verbose == 0_IK) CFC%Count%verbose = sum(CFC%Weight(1:CFC%Count%compact))
        RefinedChain%Count(RefinedChain%numRefinement)%verbose = CFC%Count%verbose

        if (present(sampleRefinementMethod)) then
            sampleRefinementMethodLowerCase = getLowerCase(sampleRefinementMethod)
            if (index(sampleRefinementMethodLowerCase,getLowerCase(BATCH_MEANS_METHOD_NAME))>0) then
                Method%isBatchMeans = .true.
            elseif (index(sampleRefinementMethodLowerCase,getLowerCase(CUTOFF_AUTOCORR_METHOD_NAME))>0 .or. index(sampleRefinementMethodLowerCase,"cutoff")>0) then
                Method%isCutoffAutoCorr = .true.
            elseif (index(sampleRefinementMethodLowerCase,getLowerCase(MAX_CUMSUM_AUTOCORR_METHOD_NAME))>0 .or. index(sampleRefinementMethodLowerCase,"cumsum")>0) then
                Method%isMaxCumSumAutoCorr = .true.
            else
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Unknown unsupported IAC computation method name: " // sampleRefinementMethod
                return
            end if
            Method%isViaCompactChain = index(sampleRefinementMethodLowerCase,"compact") > 0
            Method%isViaVerboseChain = index(sampleRefinementMethodLowerCase,"verbose") > 0
            if (.not.(Method%isViaCompactChain .or. Method%isViaVerboseChain)) then
                Method%isViaCompactChain = .true.
                Method%isViaVerboseChain = .true.
            end if
        else
            Method%isBatchMeans = .true.     ! this is the default method
        end if


        ! assign the column headers

        if (allocated(CFC%ColHeader)) then
            if (allocated(RefinedChain%ColHeader)) deallocate(RefinedChain%ColHeader)
            allocate(RefinedChain%ColHeader(0:RefinedChain%ndim))
            do i = 0, RefinedChain%ndim
                RefinedChain%ColHeader(i)%record = CFC%ColHeader(i+NUM_DEF_COL)%record
            end do
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Internal error occurred. CFC%ColHeader is NULL."
            return
        end if

        if (present(burninLoc)) then
            burninLocDefault = burninLoc
        else
            burninLocDefault = CFC%BurninLoc(CFC%Count%compact)
        end if

        RefinedChain%Count(RefinedChain%numRefinement)%compact = CFC%Count%compact - burninLocDefault + 1
        if (allocated(RefinedChain%LogFuncState)) deallocate(RefinedChain%LogFuncState)
        allocate(RefinedChain%LogFuncState(0:RefinedChain%ndim,RefinedChain%Count(RefinedChain%numRefinement)%compact))
        RefinedChain%LogFuncState(0                  ,1:RefinedChain%Count(RefinedChain%numRefinement)%compact) = CFC%LogFunc(burninLocDefault:CFC%Count%compact)
        RefinedChain%LogFuncState(1:RefinedChain%ndim,1:RefinedChain%Count(RefinedChain%numRefinement)%compact) = CFC%State(1:RefinedChain%ndim, burninLocDefault:CFC%Count%compact)

        ! check if there are more than 1 sample points in the burnin-subtracted CFC

        if (RefinedChain%Count(RefinedChain%numRefinement)%compact==0_IK) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": The size of the refined sample is zero."
            return
        elseif (RefinedChain%Count(RefinedChain%numRefinement)%compact==1_IK) then
            if (allocated(RefinedChain%Weight)) deallocate(RefinedChain%Weight)
            allocate(RefinedChain%Weight(RefinedChain%Count(RefinedChain%numRefinement)%compact))
            if (refinedChainSizeIsPresent) then
                RefinedChain%Weight = refinedChainSize
            else
                RefinedChain%Weight = 1
            end if
            if (allocated(RefinedChain%IAC)) deallocate(RefinedChain%IAC)
            if (allocated(RefinedChain%Count)) deallocate(RefinedChain%Count)
            allocate(RefinedChain%IAC(0:RefinedChain%ndim,0:0))
            allocate(RefinedChain%Count(0:0))
            RefinedChain%IAC = 0._RK
            RefinedChain%Count(RefinedChain%numRefinement)%verbose = sum(RefinedChain%Weight)
            return
        end if

        RefinedChain%Weight = CFC%Weight(burninLocDefault:CFC%Count%compact) ! Weight is intentionally separately assigned from State here

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        loopRefinement: do

            countCompact = RefinedChain%Count(RefinedChain%numRefinement)%compact

            ! set the sample weight

            if (Method%isViaCompactChain) then
                if (allocated(SampleWeight)) deallocate(SampleWeight)
            elseif (Method%isViaVerboseChain) then
                SampleWeight = RefinedChain%Weight
            endif

            ! obtain the IAC for each individual variable

            do i = 0, RefinedChain%ndim
                DumVec = RefinedChain%LogFuncState(i,1:RefinedChain%Count(RefinedChain%numRefinement)%compact)
                if (Method%isBatchMeans) then
                    RefinedChain%IAC(i,RefinedChain%numRefinement) = getBatchMeansIAC   ( np        = countCompact  &
                                                                                        , Point     = DumVec        &
                                                                                        , Weight    = SampleWeight  &
                                                                                        )
                elseif (Method%isCutoffAutoCorr) then
                    RefinedChain%IAC(i,RefinedChain%numRefinement) = getCumSumIAC       ( np        = countCompact  &
                                                                                        , Point     = DumVec        &
                                                                                        , Weight    = SampleWeight  &
                                                                                        )
                elseif (Method%isMaxCumSumAutoCorr) then
                    RefinedChain%IAC(i,RefinedChain%numRefinement) = getMaxCumSumIAC    ( np        = countCompact  &
                                                                                        , Point     = DumVec        &
                                                                                        , Weight    = SampleWeight  &
                                                                                        )
                end if
            end do

            if (refinedChainSizeIsPresent) then
                integratedAutoCorrTime = real(sum(RefinedChain%Weight),kind=RK) / real(refinedChainSize,kind=RK)
            else
                integratedAutoCorrTime = maxval( RefinedChain%IAC(0:RefinedChain%ndim,RefinedChain%numRefinement) )
            end if

            ! so far, we have computed the max IAC of the sample, but no refinement. Refine the sample only if needed.

            maxRefinementCountIsReached = RefinedChain%numRefinement==maxRefinementCount
            if (integratedAutoCorrTime<2._RK .or. maxRefinementCountIsReached) then
                if (Method%isViaCompactChain .and. Method%isViaVerboseChain) then
                    !if (maxRefinementCountIsReached) maxRefinementCount = maxRefinementCount * 2_IK
                    maxRefinementCount = maxRefinementCount * 2_IK
                    Method%isViaCompactChain = .false.
                    cycle loopRefinement
                end if
                exit loopRefinement
            end if

            ! generate the refined sample, dump it in CFC, then put it back into RefinedChain to start over again

            RefinedChain%numRefinement = RefinedChain%numRefinement + 1_IK

            ! reallocate to bigger array if nedded

            if (RefinedChain%numRefinement>maxRefinementCurrent) then
                call move_alloc( from = RefinedChain%IAC, to = DumIAC )
                call move_alloc( from = RefinedChain%Count, to = DumCount )
                maxRefinementCurrent = min( maxRefinementCurrent*2 , maxRefinementCount )
                allocate( RefinedChain%IAC( 0:RefinedChain%ndim , 0:maxRefinementCurrent ) )
                allocate( RefinedChain%Count( 0:maxRefinementCurrent ) )
                RefinedChain%IAC(0:RefinedChain%ndim,0:RefinedChain%numRefinement-1) = DumIAC
                RefinedChain%Count(0:RefinedChain%numRefinement-1) = DumCount
            end if

            if (integratedAutoCorrTime<2._RK) cycle loopRefinement ! no need for refinement. should happen only when transitioning from compact to verbose

            call refineWeightedSample   ( nd = RefinedChain%ndim                                        &
                                        , np = countCompact                                             &
                                        , skip = integratedAutoCorrTime                                 &
                                        , Sample = RefinedChain%LogFuncState                            &
                                        , Weight = RefinedChain%Weight                                  &
                                        , RefinedChain = DumCFC%State                                   &
                                        , RefinedWeight = DumCFC%Weight                                 &
                                        , PointCount = RefinedChain%Count(RefinedChain%numRefinement)   &
                                        , refinedChainSize = refinedChainSize                           &
                                        )
            RefinedChain%Weight        = DumCFC%Weight
            RefinedChain%LogFuncState  = DumCFC%State

        end do loopRefinement

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end subroutine getRefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! examples:
    !          skip: 1
    !        Weight: 5, 0, 1, 3, 1
    ! RefinedWeight: 5, 0, 1, 3, 1
    !          skip: 2
    !        Weight: 5, 0, 1, 3, 1
    ! RefinedWeight: 3, 0, 0, 2, 0
    !          skip: 3
    !        Weight: 5, 0, 1, 3, 1
    ! RefinedWeight: 2, 0, 0, 1, 1
    pure function getRefinedWeight(np,Weight,skip,refinedChainSize) result(RefinedWeight)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefinedWeight
#endif
        use Math_mod, only: getCumSum
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(in)                :: skip
        integer(IK), intent(in)             :: np, Weight(np)
        integer(IK) , intent(in) , optional :: refinedChainSize
        integer(IK)                         :: RefinedWeight(np)
        integer(IK)                         :: ip, refinedChainSizeCounter, offset
        real(RK)                            :: sumSkips, CumSumWeight(np)
        logical                             :: refinedChainSizeIsPresent
        refinedChainSizeIsPresent = present(refinedChainSize)
        if (refinedChainSizeIsPresent) refinedChainSizeCounter = 0_IK
        CumSumWeight = real(getCumSum(np,Weight),kind=RK)
        sumSkips = skip
        offset = 1_IK
        ip = offset
        RefinedWeight = 0_IK
        loopOverAllSample: do
            loopOverCurrentSample: do
                if (sumSkips>CumSumWeight(ip)) then
                    exit loopOverCurrentSample
                elseif (refinedChainSizeIsPresent) then
                    if (refinedChainSizeCounter==refinedChainSize) exit loopOverAllSample
                    refinedChainSizeCounter = refinedChainSizeCounter + 1_IK
                end if
                RefinedWeight(ip) = RefinedWeight(ip) + 1_IK
                sumSkips = sumSkips + skip
                cycle loopOverCurrentSample
            end do loopOverCurrentSample
            if (ip==np) then
                if (refinedChainSizeIsPresent) then
                    if (refinedChainSizeCounter<refinedChainSize) then
                        offset = offset + 1_IK
                        if (offset==np) offset = 1_IK
                        ip = offset
                        sumSkips = skip
                        if (offset/=1_IK) sumSkips = sumSkips + CumSumWeight(ip-1)
                        cycle loopOverAllSample
                    end if
                end if
                exit loopOverAllSample
            end if
            ip = ip + 1_IK
        end do loopOverAllSample
    end function getRefinedWeight

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine refineWeightedSample(nd,np,skip,Sample,Weight,RefinedChain,RefinedWeight,PointCount,refinedChainSize)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: refineWeightedSample
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)     , intent(in)                :: nd, np, Weight(np)   !, skip
        real(RK)        , intent(in)                :: Sample(0:nd,np), skip
        real(RK)        , intent(out), allocatable  :: RefinedChain(:,:)
        integer(IK)     , intent(out), allocatable  :: RefinedWeight(:)
        integer(IK)     , intent(in) , optional     :: refinedChainSize
        type(Count_type), intent(out)               :: PointCount
        integer(IK)                                 :: ip, ipRefined, npRefined, UpdatedWeight(np)
        UpdatedWeight = getRefinedWeight(np,Weight,skip,refinedChainSize)
        npRefined = count(UpdatedWeight>0)
        allocate( RefinedChain(0:nd,npRefined) , RefinedWeight(npRefined) )
        ipRefined = 0
        PointCount%verbose = 0
        do ip = 1, np
            if (UpdatedWeight(ip)>0) then
                ipRefined = ipRefined + 1
                RefinedChain(0:nd,ipRefined) = Sample(0:nd,ip)
                RefinedWeight(ipRefined) = UpdatedWeight(ip)
                PointCount%verbose = PointCount%verbose + RefinedWeight(ipRefined)
            end if
        end do
        PointCount%compact = npRefined
    end subroutine refineWeightedSample

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getSkip4NewSampleSize(oldSampleSize,newSampleSize) result(skip4NewSampleSize)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSkip4NewSampleSize
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: oldSampleSize,newSampleSize
        integer(IK)                 :: skip4NewSampleSize, addition, quotient
        addition = 1
        quotient = oldSampleSize / newSampleSize
        if (mod(oldSampleSize,newSampleSize)==0) addition = 0
        skip4NewSampleSize = quotient + addition
    end function getSkip4NewSampleSize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeRefinedChain(RefinedChain,sampleFileUnit,sampleFileHeaderFormat,sampleFileContentsFormat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: writeRefinedChain
#endif
        implicit none
        class(RefinedChain_type), intent(in)    :: RefinedChain
        integer(IK) , intent(in)                :: sampleFileUnit
        character(*), intent(in)                :: sampleFileHeaderFormat, sampleFileContentsFormat
        integer(IK)                             :: isample, iweight, i
        write(sampleFileUnit, sampleFileHeaderFormat) (trim(adjustl(RefinedChain%ColHeader(i)%record)),i=0,RefinedChain%ndim)
        do isample = 1, RefinedChain%Count(RefinedChain%numRefinement)%compact
            do iweight = 1, RefinedChain%Weight(isample)
                write(sampleFileUnit,sampleFileContentsFormat) RefinedChain%LogFuncState(0:RefinedChain%ndim,isample)
            end do
        end do
        flush(sampleFileUnit)
    end subroutine writeRefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function readRefinedChain(sampleFilePath,delimiter,ndim) result(RefinedChain)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: readRefinedChain
#endif
        use FileContents_mod, only: getNumRecordInFile
        use String_mod, only: String_type, num2str
        implicit none
        integer(IK) , intent(in)            :: ndim
        character(*), intent(in)            :: sampleFilePath, delimiter
        type(RefinedChain_type)             :: RefinedChain
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@readRefinedChain()"
        integer(IK)                         :: sampleFileUnit, isample, i
        type(String_type)                   :: Record

        if (allocated(Record%value)) deallocate(Record%value)
        allocate( character(99999) :: Record%value )

        RefinedChain%numRefinement = 0_IK
        RefinedChain%ndim = ndim
        allocate(RefinedChain%Count(RefinedChain%numRefinement:RefinedChain%numRefinement))

        ! find the number of lines in the sample file

        call getNumRecordInFile(filePath=sampleFilePath,numRecord=RefinedChain%Count(RefinedChain%numRefinement)%verbose,Err=RefinedChain%Err)
        if (RefinedChain%Err%occurred) then
            RefinedChain%Err%msg = PROCEDURE_NAME // RefinedChain%Err%msg
            return
        end if

        RefinedChain%Count(RefinedChain%numRefinement)%verbose = RefinedChain%Count(RefinedChain%numRefinement)%verbose - 1_IK ! remove header from the count
        allocate( RefinedChain%LogFuncState(RefinedChain%Count(RefinedChain%numRefinement)%verbose,0:ndim) )

        open( newunit = sampleFileUnit &
            , file = sampleFilePath &
            , status = "old" &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
            , SHARED &
#endif
            )

        ! read header

        read(sampleFileUnit,"(A)") Record%value
        Record%Parts = Record%SplitStr(Record%value,delimiter,Record%nPart)
        if (Record%nPart/=ndim+1_IK) then
            RefinedChain%Err%occurred = .true.
            RefinedChain%Err%msg = PROCEDURE_NAME // ": The number of column headers ("//num2str(Record%nPart)//") is not equal to ndim + 1: "//num2str(ndim+1_IK)
            return
        end if

        allocate(RefinedChain%ColHeader(0:ndim))
        do i = 0, ndim
            RefinedChain%ColHeader(i)%record = Record%Parts(i+1)%record
        end do

        ! read contents

        do isample = 1, RefinedChain%Count(RefinedChain%numRefinement)%verbose
            read(sampleFileUnit, "(A)") Record%value
            Record%Parts = Record%SplitStr(trim(adjustl(Record%value)),delimiter)
            do i = 0, ndim
                read(Record%Parts(i+1)%record,*) RefinedChain%LogFuncState(isample,i)
            end do
        end do

        close(sampleFileUnit)

    end function readRefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module ParaMCMCRefinedChain_mod