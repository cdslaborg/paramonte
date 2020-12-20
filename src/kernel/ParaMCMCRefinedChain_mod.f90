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

!>  \brief This module contains the classes and procedures for refining MCMC output chains.
!>  \author Amir Shahmoradi

module ParaMCMCRefinedChain_mod

    use SpecMCMC_SampleRefinementMethod_mod, only: BATCH_MEANS_METHOD_NAME
    use SpecMCMC_SampleRefinementMethod_mod, only: CUTOFF_AUTOCORR_METHOD_NAME
    use SpecMCMC_SampleRefinementMethod_mod, only: MAX_CUMSUM_AUTOCORR_METHOD_NAME
    use ParaMonteChainFileContents_mod, only: Count_type
    use JaggedArray_mod, only: CharVec_type
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "\paramCMCRefinedChain_mod"

    !> The `RefinedChain_type` class.
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
       logical :: isMax = .false.
       logical :: isMed = .false.
       logical :: isMin = .false.
       logical :: isAvg = .false.
       logical :: isBatchMeans = .false.
       logical :: isCutoffAutoCorr = .false.
       logical :: isViaCompactChain = .true.
       logical :: isViaVerboseChain = .true.
       logical :: isMaxCumSumAutoCorr = .false.
   end type Method_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the refined Markov chain, given the input Markov chain and its specifications.
    !> This procedure is a method of the [RefinedChain_type](@ref refinedchain_type) class.
    !>
    !> \param[inout] RefinedChain           :   An object of class [RefinedChain_type](@ref refinedchain_type).
    !> \param[inout] CFC                    :   An object of type [ChainFileContents_type](@ref paramontechainfilecontents_mod::chainfilecontents_type)
    !!                                          containing the Markov chain.
    !> \param[out]   Err                    :   An object of class [Err_type](@ref err_mod::err_type) indicating whether any error has occurred or not.
    !> \param[in]    burninLoc              :   The estimated location of burnin point in the Markov chain (**optional**).
    !!                                          If not provided, it will be extracted from the components of the input `CFC`.
    !> \param[in]    refinedChainSize       :   The requested refined sample size (**optional**). If the size of the refined sample is given as input,
    !!                                          then the requested sample is directly generated based on the input size.
    !> \param[in]    sampleRefinementCount  :   The maximum number of times the sample can be refined (**optional**, default = `Infinity`).
    !!                                      :   For example, if set to 1, then only one round of refinement will be performed on the Markov chain.
    !> \param[in]    sampleRefinementMethod :   The requested method of refining the sample (**optional**, default = "BatchMeans").
    subroutine getRefinedChain  ( RefinedChain              &
                                , CFC                       &
                                , Err                       &
                                , burninLoc                 &
                                , refinedChainSize          &
                                , sampleRefinementCount     &
                                , sampleRefinementMethod    &
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefinedChain
#endif

        use, intrinsic :: iso_fortran_env, only: output_unit
        use ParaMonteChainFileContents_mod, only: ChainFileContents_type, NUM_DEF_COL
        use CrossCorr_mod, only: getBatchMeansIAC, getCumSumIAC, getMaxCumSumIAC
        use String_mod, only: getLowerCase, replaceStr
        use Sort_mod, only: getMedian

        implicit none

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@getRefinedChain()"

        class(RefinedChain_type)    , intent(inout)             :: RefinedChain
        type(ChainFileContents_type), intent(inout)             :: CFC
        type(Err_type)              , intent(out)               :: Err
        character(*)                , intent(in), optional      :: sampleRefinementMethod
        integer(IK)                 , intent(in), optional      :: burninLoc, refinedChainSize, sampleRefinementCount
        integer(IK)                                             :: burninLocDefault, i, countCompact, ndimPlusOne
        integer(IK)                                             :: maxRefinementCurrent, maxRefinementCount
        real(RK)                                                :: integratedAutoCorrTime, ndimPlusOneInverse
        logical                                                 :: refinedChainSizeIsPresent
        logical                                                 :: maxRefinementCountIsReached
        character(:)    , allocatable                           :: sampleRefinementMethodDefault, sampleRefinementMethodDefaultLowerCase, srmethod
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

        ! define the AIC computation method

        sampleRefinementMethodDefault = BATCH_MEANS_METHOD_NAME; if (present(sampleRefinementMethod)) sampleRefinementMethodDefault = sampleRefinementMethod
        sampleRefinementMethodDefaultLowerCase = getLowerCase(sampleRefinementMethodDefault)

        if (index(sampleRefinementMethodDefaultLowerCase,getLowerCase(BATCH_MEANS_METHOD_NAME))>0) then
            Method%isBatchMeans = .true.
        elseif (index(sampleRefinementMethodDefaultLowerCase,getLowerCase(CUTOFF_AUTOCORR_METHOD_NAME))>0 .or. index(sampleRefinementMethodDefaultLowerCase,"cutoff")>0) then
            Method%isCutoffAutoCorr = .true.
        elseif (index(sampleRefinementMethodDefaultLowerCase,getLowerCase(MAX_CUMSUM_AUTOCORR_METHOD_NAME))>0 .or. index(sampleRefinementMethodDefaultLowerCase,"cumsum")>0) then
            Method%isMaxCumSumAutoCorr = .true.
        else
            ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Unknown unsupported IAC computation method name: " // sampleRefinementMethodDefault
            return
            ! LCOV_EXCL_STOP
        end if

        ! define the chain types to use for the AIC computation

        Method%isViaCompactChain = index(sampleRefinementMethodDefaultLowerCase,"compact") > 0
        Method%isViaVerboseChain = index(sampleRefinementMethodDefaultLowerCase,"verbose") > 0
        if (.not.(Method%isViaCompactChain .or. Method%isViaVerboseChain)) then
            Method%isViaCompactChain = .true.
            Method%isViaVerboseChain = .true.
        end if

        ! define the statistic to use for the AIC computation

        srmethod = replaceStr(string=sampleRefinementMethodDefaultLowerCase,search=" ",substitute="-")
        Method%isAvg =  index(srmethod,"-avg") > 0 .or. index(srmethod,"-average") > 0
        Method%isMed =  index(srmethod,"-med") > 0 .or. index(srmethod,"-median") > 0
        Method%isMin =  index(srmethod,"-min") > 0 .or. index(srmethod,"-minimum") > 0
        Method%isMax =  index(srmethod,"-max") > 0 .or. index(srmethod,"-maximum") > 0
        if ( Method%isAvg .and. (Method%isMed .or. Method%isMax .or. Method%isMin) ) Err%occurred = .true.
        if ( Method%isMed .and. (Method%isMax .or. Method%isMin) ) Err%occurred = .true.
        if ( Method%isMax .and. Method%isMin ) Err%occurred = .true.
        if (.not.(Method%isAvg .or. Method%isMed .or. Method%isMin .or. Method%isMax)) Method%isAvg = .true. ! default method of AIC summarization.

        if (Method%isAvg) ndimPlusOneInverse = 1._RK / (RefinedChain%ndim + 1_IK)
        if (Method%isMed) ndimPlusOne = RefinedChain%ndim + 1_IK

        ! assign the column headers

        if (allocated(CFC%ColHeader)) then
            if (allocated(RefinedChain%ColHeader)) deallocate(RefinedChain%ColHeader)
            allocate(RefinedChain%ColHeader(0:RefinedChain%ndim))
            do i = 0, RefinedChain%ndim
                RefinedChain%ColHeader(i)%record = CFC%ColHeader(i+NUM_DEF_COL)%record
            end do
        else
            ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Internal error occurred. CFC%ColHeader is NULL."
            return
            ! LCOV_EXCL_STOP
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
            ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": The size of the refined sample is zero."
            return
            ! LCOV_EXCL_STOP
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
                                                                                        , Point     = DumVec        & ! LCOV_EXCL_LINE
                                                                                        , Weight    = SampleWeight  & ! LCOV_EXCL_LINE
                                                                                        )
                elseif (Method%isCutoffAutoCorr) then
                    RefinedChain%IAC(i,RefinedChain%numRefinement) = getCumSumIAC       ( np        = countCompact  &
                                                                                        , Point     = DumVec        & ! LCOV_EXCL_LINE
                                                                                        , Weight    = SampleWeight  & ! LCOV_EXCL_LINE
                                                                                        )
                elseif (Method%isMaxCumSumAutoCorr) then
                    RefinedChain%IAC(i,RefinedChain%numRefinement) = getMaxCumSumIAC    ( np        = countCompact  &
                                                                                        , Point     = DumVec        & ! LCOV_EXCL_LINE
                                                                                        , Weight    = SampleWeight  & ! LCOV_EXCL_LINE
                                                                                        )
                end if
            end do

            if (refinedChainSizeIsPresent) then
                integratedAutoCorrTime = real(sum(RefinedChain%Weight),kind=RK) / real(refinedChainSize,kind=RK)
            else
                if (Method%isAvg) then
                    integratedAutoCorrTime = sum( RefinedChain%IAC(0:RefinedChain%ndim,RefinedChain%numRefinement) ) * ndimPlusOneInverse
                elseif (Method%isMed) then
                    call getMedian(lenArray=ndimPlusOne,Array=RefinedChain%IAC(0:RefinedChain%ndim,RefinedChain%numRefinement),median=integratedAutoCorrTime,Err=Err)
                    if (Err%occurred) then
                        ! LCOV_EXCL_START
                        Err%msg = PROCEDURE_NAME//Err%msg
                        return
                        ! LCOV_EXCL_STOP
                    end if
                elseif (Method%isMax) then
                    integratedAutoCorrTime = maxval( RefinedChain%IAC(0:RefinedChain%ndim,RefinedChain%numRefinement) )
                elseif (Method%isMin) then
                    integratedAutoCorrTime = minval( RefinedChain%IAC(0:RefinedChain%ndim,RefinedChain%numRefinement) )
                end if
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

            call refineWeightedSample   ( nd = RefinedChain%ndim                                        & ! LCOV_EXCL_LINE
                                        , np = countCompact                                             & ! LCOV_EXCL_LINE
                                        , skip = integratedAutoCorrTime                                 & ! LCOV_EXCL_LINE
                                        , Sample = RefinedChain%LogFuncState                            & ! LCOV_EXCL_LINE
                                        , Weight = RefinedChain%Weight                                  & ! LCOV_EXCL_LINE
                                        , RefinedChain = DumCFC%State                                   & ! LCOV_EXCL_LINE
                                        , RefinedWeight = DumCFC%Weight                                 & ! LCOV_EXCL_LINE
                                        , PointCount = RefinedChain%Count(RefinedChain%numRefinement)   &
                                        , refinedChainSize = refinedChainSize                           & ! LCOV_EXCL_LINE
                                        )
            RefinedChain%Weight        = DumCFC%Weight
            RefinedChain%LogFuncState  = DumCFC%State

        end do loopRefinement

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end subroutine getRefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the refined vector of weights of the vector of weights of a weighted Markov chain.
    !>
    !> \param[in]   np                  :   The number of elements of the `Weight` vector.
    !> \param[in]   Weight              :   The input vector of weights.
    !> \param[in]   skip                :   The size of the jumps that have to be made through the weighted Markov chain.
    !> \param[in]   refinedChainSize    :   The requested refined sample size (**optional**). If present, then the refined chain (represented by the
    !>                                  :   vector `Weight`) will be refined such that the resulting refined chain has the size `refinedChainSize`.
    !>
    !> \return
    !> `RefinedWeight` : An array of size `np`, whose elements indicate which points are present in the final refined chain.\n
    !> Examples:
    !> ```
    !>          skip: 1
    !>        Weight: 5, 0, 1, 3, 1
    !> RefinedWeight: 5, 0, 1, 3, 1
    !>          skip: 2
    !>        Weight: 5, 0, 1, 3, 1
    !> RefinedWeight: 3, 0, 0, 2, 0
    !>          skip: 3
    !>        Weight: 5, 0, 1, 3, 1
    !> RefinedWeight: 2, 0, 0, 1, 1
    !> ```
    pure function getRefinedWeight(np,Weight,skip,refinedChainSize) result(RefinedWeight)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
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
                        cycle loopOverAllSample ! LCOV_EXCL_LINE
                    end if
                end if
                exit loopOverAllSample
            end if
            ip = ip + 1_IK
        end do loopOverAllSample
    end function getRefinedWeight

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Refined an input weighted sample according to the new requested weights.
    !>
    !> \param[in]   nd                  :   The number of dimensions of the input `Sample(0:nd,np)`.
    !> \param[in]   np                  :   The number of points in the input `Sample(0:nd,np)`.
    !> \param[in]   skip                :   The jump size with which the input chain has to be refined.
    !> \param[in]   Sample              :   The input 2-dimensional array of sampled states which has to be refined.
    !> \param[in]   Weight              :   The weights of the sampled points.
    !> \param[out]  RefinedChain        :   The refined array.
    !> \param[out]  RefinedWeight       :   The vector of refined weights corresponding to the output refined array.
    !> \param[out]  PointCount          :   An object of derived type [Count_type](@ref paramontechainfilecontents_mod::count_type)
    !>                                      containing the number of points in the refined sample.
    !> \param[in]   refinedChainSize    :   The requested refined sample size (**optional**). If the size of the refined sample is given as input,
    !>                                      then the requested sample is directly generated based on the input size.
    pure subroutine refineWeightedSample(nd,np,skip,Sample,Weight,RefinedChain,RefinedWeight,PointCount,refinedChainSize)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
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
        integer(IK)                                 :: ip, ipRefined, npRefined, UpdatedWeight(np) ! LCOV_EXCL_LINE
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

    !> Return the best skip size through a Markov chain to refined it to the optimal requested size.
    !>
    !> \param[in]   oldSampleSize   :   The original size of the Markov chain.
    !> \param[in]   newSampleSize   :   The final desired size of the refined sample.
    !>
    !> \return
    !> `skip4NewSampleSize` : The computed skip size.
    !>
    !> \warning
    !> The condition `oldSampleSize >= newSampleSize` must always hold,
    !> otherwise a negative value for `skip4NewSampleSize` will be returned to indicate the occurrence of an error.
    pure function getSkip4NewSampleSize(oldSampleSize,newSampleSize) result(skip4NewSampleSize)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSkip4NewSampleSize
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: oldSampleSize,newSampleSize
        integer(IK)                 :: skip4NewSampleSize, addition, quotient
        if (oldSampleSize < newSampleSize) then
            skip4NewSampleSize = -1_IK
            return
        end if
        addition = 1
        quotient = oldSampleSize / newSampleSize
        if (mod(oldSampleSize,newSampleSize)==0) addition = 0
        skip4NewSampleSize = quotient + addition
    end function getSkip4NewSampleSize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Write the computed refined chain to the specified output file.
    !>
    !> \param[in]   RefinedChain                :   An object of class [RefinedChain_type](@ref refinedchain_type)
    !>                                              containing the refined sample to be written to the output file.
    !> \param[in]   sampleFileUnit              :   The unit of the file to which the sample must be written.
    !> \param[in]   sampleFileHeaderFormat      :   The IO format of the header of the sample file.
    !> \param[in]   sampleFileContentsFormat    :   The IO format of the contents (sampled states) in the sample file.
    subroutine writeRefinedChain(RefinedChain,sampleFileUnit,sampleFileHeaderFormat,sampleFileContentsFormat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
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
    end subroutine writeRefinedChain ! LCOV_EXCL_LINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Write the computed refined chain to the specified output file.
    !>
    !> \param[in]   sampleFilePath      :   The path to the input chain file that must be read.
    !> \param[in]   delimiter           :   The delimiter used in the file.
    !> \param[in]   ndim                :   The number of dimensions of the sampled states in the sample file.
    !>                                      This is basically, the size of the domain of the objective function.
    !> \param[in]   tenabled            :   An optional input logical value standing for `transpose-enabled`. If `.false.`,
    !>                                      the input data will be naturally read according to Fortran column-wise data storage
    !>                                      rule as a matrix of rank `0:nd * 1:np`. If `.false.`, the input sample file will be
    !>                                      read as a matrix of rank `1:np * 0:nd`. Note that `np` represents the number of rows
    !>                                      in the files (that is, the number of sampled points, whereas `nd` represents the
    !>                                      number of columns in the input file (**optional**, default = `.false.`).
    !>
    !> \return
    !> `RefinedChain` : An object of class [RefinedChain_type](@ref refinedchain_type) containing
    !>                  the sampled states read from the specified input file.
    function readRefinedChain(sampleFilePath,delimiter,ndim,tenabled) result(RefinedChain)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: readRefinedChain
#endif
        use FileContents_mod, only: getNumRecordInFile
        use String_mod, only: String_type, num2str
        implicit none
        integer(IK) , intent(in)            :: ndim
        character(*), intent(in)            :: sampleFilePath, delimiter
        logical     , intent(in), optional  :: tenabled
        type(RefinedChain_type)             :: RefinedChain
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@readRefinedChain()"
        integer(IK)                         :: sampleFileUnit, isample, i
        logical                             :: tenabledDefault
        type(String_type)                   :: Record

        if (allocated(Record%value)) deallocate(Record%value) ! LCOV_EXCL_LINE
        allocate( character(99999) :: Record%value )

        RefinedChain%ndim = ndim
        RefinedChain%numRefinement = 0_IK
        allocate(RefinedChain%Count(RefinedChain%numRefinement:RefinedChain%numRefinement)) ! allocate just the zeroth level of `RefinedChain`.

        ! find the number of lines in the sample file

        call getNumRecordInFile(filePath=sampleFilePath, numRecord=RefinedChain%Count(RefinedChain%numRefinement)%verbose, Err=RefinedChain%Err)
        if (RefinedChain%Err%occurred) then
            ! LCOV_EXCL_START
            RefinedChain%Err%msg = PROCEDURE_NAME // RefinedChain%Err%msg
            return
            ! LCOV_EXCL_STOP
        end if

        RefinedChain%Count(RefinedChain%numRefinement)%verbose = RefinedChain%Count(RefinedChain%numRefinement)%verbose - 1_IK ! remove header from the count

        open( newunit = sampleFileUnit &
            , file = sampleFilePath &
            , status = "old" &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED &
#endif
            )

        ! read header

        read(sampleFileUnit,"(A)") Record%value
        Record%Parts = Record%split(Record%value,delimiter,Record%nPart)
        if (Record%nPart/=ndim+1_IK) then
            ! LCOV_EXCL_START
            RefinedChain%Err%occurred = .true.
            RefinedChain%Err%msg = PROCEDURE_NAME // ": The number of header columns ("//num2str(Record%nPart)//") is not equal to ndim + 1: "//num2str(ndim+1_IK)
            return
            ! LCOV_EXCL_STOP
        end if

        allocate(RefinedChain%ColHeader(0:ndim))
        do i = 0, ndim
            RefinedChain%ColHeader(i)%record = Record%Parts(i+1)%record
        end do

        ! read contents

        tenabledDefault = .false.
        if (present(tenabled)) tenabledDefault = tenabled

        if (tenabledDefault) then

            allocate( RefinedChain%LogFuncState(RefinedChain%Count(RefinedChain%numRefinement)%verbose, 0:ndim) )
            do isample = 1, RefinedChain%Count(RefinedChain%numRefinement)%verbose
                read(sampleFileUnit, "(A)") Record%value
                Record%Parts = Record%split(trim(adjustl(Record%value)),delimiter)
                do i = 0, ndim
                    read(Record%Parts(i+1)%record,*) RefinedChain%LogFuncState(isample,i)
                end do
            end do

        else

            allocate( RefinedChain%LogFuncState(0:ndim, RefinedChain%Count(RefinedChain%numRefinement)%verbose) )
            do isample = 1, RefinedChain%Count(RefinedChain%numRefinement)%verbose
                read(sampleFileUnit, "(A)") Record%value
                Record%Parts = Record%split(trim(adjustl(Record%value)),delimiter)
                do i = 0, ndim
                    read(Record%Parts(i+1)%record,*) RefinedChain%LogFuncState(i,isample)
                end do
            end do

        end if

        close(sampleFileUnit)

    end function readRefinedChain

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module ParaMCMCRefinedChain_mod ! LCOV_EXCL_LINE