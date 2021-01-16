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

!> \brief
!> This file implements the body of the Proposal modules of the `ParaDRAM` and `ParaDISE` samplers.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

#if defined UNIFORM

#define GET_RANDOM_PROPOSAL getRandMVU

#elif defined NORMAL

#define GET_RANDOM_PROPOSAL getRandMVN

#else

#error "Unknown Proposal model in ParaDXXX_Proposal_mod.inc.f90"

#endif


#if defined PARADRAM

#define ParaXXXX ParaDRAM
    use ParaDRAM_ProposalAbstract_mod, only: ProposalAbstract_type

#elif defined PARADISE

#define ParaXXXX ParaDISE
    use ParaDISE_ProposalAbstract_mod, only: ProposalAbstract_type

#endif

    use ParaMonte_mod, only: Image_type
    use Constants_mod, only: IK, RK, PMSM
    use String_mod, only: IntStr_type
    use Err_mod, only: Err_type

    implicit none

    private
    public :: Proposal_type

#if defined NORMAL
    character(*), parameter         :: MODULE_NAME = "@"//PMSM%ParaXXXX//"ProposalNormal_mod"
#elif defined UNIFORM
    character(*), parameter         :: MODULE_NAME = "@"//PMSM%ParaXXXX//"ProposalUniform_mod"
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type                            :: AccRate_type
        real(RK)                    :: sumUpToLastUpdate
        real(RK)                    :: target
    end type AccRate_type

    !> The `Proposal_type` class.
    type, extends(ProposalAbstract_type) :: Proposal_type
       integer(IK)                  :: logFileUnit
       integer(IK)                  :: restartFileUnit
       ! the following are made components for the sake of thread-safe save attribute and the restart file generation.
       integer(IK)                  :: sampleSizeOld
       real(RK)                     :: logSqrtDetOld
       real(RK)                     :: adaptiveScaleFactorSq
       real(RK)     , allocatable   :: MeanOld(:)
       !type(AccRate_type)          :: AccRate
    contains
        procedure   , pass          :: getNew
#if defined PARADISE
        procedure   , nopass        :: getLogProb
#endif
        procedure   , pass          :: doAdaptation
       !procedure   , nopass        :: readRestartFile
       !procedure   , nopass        :: writeRestartFile
#if defined CAF_ENABLED || defined MPI_ENABLED
        procedure   , nopass        :: bcastAdaptation
#endif
    end type Proposal_type

    interface Proposal_type
        module procedure :: constructProposalSymmetric
    end interface Proposal_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Covariance Matrix of the proposal distribution. Last index belongs to delayed rejection

#if defined CAF_ENABLED
    real(RK)        , save  , allocatable   :: comv_CholDiagLower(:,:,:)[:]
#else
    real(RK)        , save  , allocatable   :: comv_CholDiagLower(:,:,:)
#endif
    real(RK)        , save  , allocatable   :: mv_logSqrtDetInvCovMat(:)
    real(RK)        , save  , allocatable   :: mv_InvCovMat(:,:,:)

#if defined MPI_ENABLED
    integer(IK)     , save                  :: mc_ndimSqPlusNdim
#endif
    type(Image_type), save                  :: mc_Image
    integer(IK)     , save                  :: mc_ndim
    logical         , save                  :: mc_scalingRequested
    real(RK)        , save                  :: mc_defaultScaleFactorSq
    integer(IK)     , save                  :: mc_DelayedRejectionCount
    integer(IK)     , save                  :: mc_MaxNumDomainCheckToWarn
    integer(IK)     , save                  :: mc_MaxNumDomainCheckToStop
    logical         , save                  :: mc_delayedRejectionRequested
    real(RK)        , save                  :: mc_ndimInverse
    real(RK)        , save                  :: mc_targetAcceptanceRate
    real(RK)        , save                  :: mc_TargetAcceptanceRateLimit(2)
   !real(RK)        , save                  :: mc_maxScaleFactor != 2._RK
   !real(RK)        , save                  :: mc_maxScaleFactorSq != mc_maxScaleFactor**2
    real(RK)        , save  , allocatable   :: mc_DelayedRejectionScaleFactorVec(:)
    real(RK)        , save  , allocatable   :: mc_DomainLowerLimitVec(:)
    real(RK)        , save  , allocatable   :: mc_DomainUpperLimitVec(:)
    logical         , save  , allocatable   :: mc_isAsciiRestartFileFormat
    logical         , save  , allocatable   :: mc_isBinaryRestartFileFormat
    character(:)    , save  , allocatable   :: mc_MaxNumDomainCheckToWarnMsg
    character(:)    , save  , allocatable   :: mc_MaxNumDomainCheckToStopMsg
    character(:)    , save  , allocatable   :: mc_negativeTotalVariationMsg
    character(:)    , save  , allocatable   :: mc_restartFileFormat
    character(:)    , save  , allocatable   :: mc_methodBrand
    character(:)    , save  , allocatable   :: mc_methodName
#if defined UNIFORM
    real(RK)        , save  , allocatable   :: mc_negLogVolUnitBall
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> This is the constructor of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> Return an object of [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes given the input simulation characteristics.
    !>
    !> \remark
    !> This interface madness is a result of the internal compiler bug in GFortran as of Jan 2020, which diagnoses a `ParaXXXX_type`
    !> argument as circular dependency due to this constructor appearing in the type-bound setup procedure of `ParaXXXX_type`.
    !> Intel does not complain. Until GFortran comes up with a fix, we have to live with this interface.
    function constructProposalSymmetric ( ndim &
                                        , SpecBase &
                                        , SpecMCMC &
                                        , SpecDRAM &
                                        , Image &
                                        , name &
                                        , brand &
                                        , LogFile &
                                        , RestartFile &
                                        , isFreshRun &
                                        ) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalSymmetric
#endif
        use Constants_mod, only: IK, RK, NULL_RK
        use ParaMonte_mod, only: Image_type
        use ParaMonte_mod, only: LogFile_type
        use ParaMonte_mod, only: RestartFile_type
        use SpecBase_mod, only: SpecBase_type
        use SpecMCMC_mod, only: SpecMCMC_type
        use SpecDRAM_mod, only: SpecDRAM_type
        use Matrix_mod, only: getCholeskyFactor
        use String_mod, only: num2str
        use Err_mod, only: abort
#if defined UNIFORM
        use Math_mod, only: getLogVolUnitBall
#endif
        implicit none

        integer(IK)             , intent(in)    :: ndim
        type(SpecBase_type)     , intent(in)    :: SpecBase
        type(SpecMCMC_type)     , intent(in)    :: SpecMCMC
        type(SpecDRAM_type)     , intent(in)    :: SpecDRAM
        type(Image_type)        , intent(in)    :: Image
        character(*)            , intent(in)    :: name
        character(*)            , intent(in)    :: brand
        type(LogFile_type)      , intent(in)    :: LogFile
        type(RestartFile_type)  , intent(in)    :: RestartFile
        logical                                 :: isFreshRun

        type(Proposal_type)                     :: self

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@constructProposalSymmetric()"
        integer                                 :: i, j

        self%Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup sampler update global save variables
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (allocated(self%MeanOld)) deallocate(self%MeanOld); allocate(self%MeanOld(ndim))
        self%MeanOld(1:ndim) = SpecMCMC%StartPointVec%Val
        self%logSqrtDetOld   = NULL_RK
        self%sampleSizeOld   = 1_IK
        self%adaptiveScaleFactorSq = 1._RK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup general proposal specifications
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MPI_ENABLED
        mc_ndimSqPlusNdim                   = ndim*(ndim+1_IK)
#endif
        mc_ndim                             = ndim
        mc_DomainLowerLimitVec              = SpecBase%DomainLowerLimitVec%Val
        mc_DomainUpperLimitVec              = SpecBase%DomainUpperLimitVec%Val
        mc_DelayedRejectionScaleFactorVec   = SpecDRAM%delayedRejectionScaleFactorVec%Val
       !mc_isNormal                         = SpecMCMC%ProposalModel%isNormal
        mc_Image                            = Image
        mc_methodName                       = name
        mc_methodBrand                      = brand
        self%logFileUnit                    = LogFile%unit
        self%restartFileUnit                = RestartFile%unit
        mc_restartFileFormat                = RestartFile%format
        mc_isBinaryRestartFileFormat        = SpecBase%RestartFileFormat%isBinary
        mc_isAsciiRestartFileFormat         = SpecBase%RestartFileFormat%isAscii
        mc_defaultScaleFactorSq             = SpecMCMC%ScaleFactor%val**2
       !self%AccRate%sumUpToLastUpdate  = 0._RK
        mc_maxNumDomainCheckToWarn          = SpecBase%MaxNumDomainCheckToWarn%val
        mc_maxNumDomainCheckToStop          = SpecBase%MaxNumDomainCheckToStop%val
        mc_delayedRejectionCount            = SpecDRAM%DelayedRejectionCount%val
        mc_delayedRejectionRequested        = mc_DelayedRejectionCount > 0_IK
        mc_scalingRequested                 = SpecBase%TargetAcceptanceRate%scalingRequested
        mc_ndimInverse                      = 1._RK/real(ndim,kind=RK)
        mc_TargetAcceptanceRateLimit        = SpecBase%TargetAcceptanceRate%Val
        mc_targetAcceptanceRate             = 0.234_RK
        if (mc_scalingRequested) mc_targetAcceptanceRate = sum(mc_TargetAcceptanceRateLimit) / 2._RK
       !mc_maxScaleFactorSq                 = 4._RK**mc_ndimInverse
       !mc_maxScaleFactor                   = sqrt(mc_maxScaleFactorSq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! setup ProposalSymmetric specifications
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! setup covariance matrix

        if (allocated(mv_InvCovMat)) deallocate(mv_InvCovMat)
        allocate( mv_InvCovMat(ndim,0:ndim,0:mc_DelayedRejectionCount) )
        if (allocated(mv_logSqrtDetInvCovMat)) deallocate(mv_logSqrtDetInvCovMat)
        allocate( mv_logSqrtDetInvCovMat(0:mc_DelayedRejectionCount) )

        if (allocated(comv_CholDiagLower)) deallocate(comv_CholDiagLower)
#if defined CAF_ENABLED
        ! on the second dimension, the zeroth index refers to the Diagonal elements of the Cholesky lower triangular matrix
        ! This rearrangement was done for more efficient communication of the matrix across processes.
        allocate( comv_CholDiagLower(ndim,0:ndim,0:mc_DelayedRejectionCount)[*] )
#else
        allocate( comv_CholDiagLower(ndim,0:ndim,0:mc_DelayedRejectionCount) )
#endif

        comv_CholDiagLower(1:ndim,1:ndim,0) = SpecMCMC%ProposalStartCovMat%Val

        ! Now scale the covariance matrix

        do j = 1, ndim
            do i = 1, j
                comv_CholDiagLower(i,j,0) = comv_CholDiagLower(i,j,0) * mc_defaultScaleFactorSq
            end do
        end do

        ! Now get the Cholesky Factor of the Covariance Matrix. Lower comv_CholDiagLower will be the CholFac

        call getCholeskyFactor( ndim, comv_CholDiagLower(:,1:ndim,0), comv_CholDiagLower(1:ndim,0,0) ) ! The `:` instead of `1:ndim` avoids temporary array creation.
        if (comv_CholDiagLower(1,0,0)<0._RK) then
        ! LCOV_EXCL_START
            self%Err%msg = mc_Image%name // PROCEDURE_NAME // ": Singular input covariance matrix by user was detected. This is strange.\nCovariance matrix lower triangle:"
            do j = 1, ndim
                do i = 1, j
                    self%Err%msg = self%Err%msg // "\n" // num2str(comv_CholDiagLower(1:i,j,0))
                end do
            end do
            self%Err%occurred = .true.
            call abort( Err = self%Err, prefix = mc_methodBrand, newline = "\n", outputUnit = self%logFileUnit )
            return
        ! LCOV_EXCL_STOP
        end if

        if (mc_delayedRejectionRequested) call updateDelRejCholDiagLower()
        call getInvCovMat()
        self%logSqrtDetOld = sum(log( comv_CholDiagLower(1:ndim,0,0) ))

        ! Scale the higher-stage delayed-rejection Cholesky Lower matrices

        if (mc_delayedRejectionRequested) call updateDelRejCholDiagLower()

        ! This will be used for Domain boundary checking during the simulation

        mc_negativeTotalVariationMsg = MODULE_NAME//"@doAdaptation(): Non-positive adaptation measure detected, possibly due to round-off error: "
        mc_MaxNumDomainCheckToWarnMsg = MODULE_NAME//"@getNew(): "//num2str(mc_MaxNumDomainCheckToWarn)//&
                                        " proposals were drawn out of the objective function's Domain without any acceptance."
        mc_MaxNumDomainCheckToStopMsg = MODULE_NAME//"@getNew(): "//num2str(mc_MaxNumDomainCheckToStop)//&
                                        " proposals were drawn out of the objective function's Domain. As per the value set for the&
                                        & simulation specification variable 'maxNumDomainCheckToStop', "//mc_methodName//" will abort now."

#if defined UNIFORM
        if (ndim>1_IK) then
            mc_negLogVolUnitBall = -getLogVolUnitBall(ndim)
        else
            mc_negLogVolUnitBall = 0._RK
        end if
#endif

        ! read/write the first entry of the restart file

        if (mc_Image%isLeader) then
            block
                real(RK) :: meanAccRateSinceStart
                if (isFreshRun) then
                    if (mc_isBinaryRestartFileFormat) then
                        call writeRestartFileBinary(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = 1._RK)
                    else
                        call writeRestartFileAscii  ( restartFileUnit = self%restartFileUnit & ! LCOV_EXCL_LINE
                                                    , meanAccRateSinceStart = 1._RK & ! LCOV_EXCL_LINE
                                                    , sampleSizeOld = self%sampleSizeOld & ! LCOV_EXCL_LINE
                                                    , logSqrtDetOld = self%logSqrtDetOld & ! LCOV_EXCL_LINE
                                                    , adaptiveScaleFactorSq = self%adaptiveScaleFactorSq & ! LCOV_EXCL_LINE
                                                    , MeanOld = self%MeanOld & ! LCOV_EXCL_LINE
                                                    )
                    end if
                else
                    if (mc_isBinaryRestartFileFormat) then
                        call readRestartFileBinary(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = meanAccRateSinceStart)
                    else
                        call readRestartFileAscii(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = meanAccRateSinceStart)
                    end if
                end if
            end block
        end if

    end function constructProposalSymmetric

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> Return a newly sampled state to the kernel routine for acceptance check.
    !>
    !>
    !> @param[in]   nd          :   The number of dimensions of the objective function.
    !> @param[in]   counterDRS  :   The Delayed Rejection Stage counter.
    !> @param[in]   StateOld    :   The old sampled state.
    !>
    !> \return
    !> `StateNew` : The newly sampled state.
    subroutine getNew   ( self          &
                        , nd            &
                        , counterDRS    &
                        , StateOld      &
                        , StateNew      &
                        )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNew
#endif
        use Statistics_mod, only: GET_RANDOM_PROPOSAL
        use Constants_mod, only: IK, RK
        use Err_mod, only: warn, abort
        use ParaMonteLogFunc_mod, only: getLogFunc_proc

        implicit none

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getNew()"

        class(Proposal_type), intent(inout) :: self
        integer(IK), intent(in)             :: nd
        integer(IK), intent(in)             :: counterDRS
        real(RK)   , intent(in)             :: StateOld(nd)
        real(RK)   , intent(out)            :: StateNew(nd)
        integer(IK)                         :: domainCheckCounter

        domainCheckCounter = 0_IK

        loopBoundaryCheck: do ! Check for the support Region consistency:
#if defined UNIFORM || defined NORMAL
            StateNew(1:nd) = GET_RANDOM_PROPOSAL( nd                                    & ! LCOV_EXCL_LINE
                                                , StateOld                              & ! LCOV_EXCL_LINE
                                                ! ATTN: The colon index in place of 1:nd below avoids the temporary array creation
                                                , comv_CholDiagLower(:,1:nd,counterDRS) & ! LCOV_EXCL_LINE
                                                , comv_CholDiagLower(1:nd,0,counterDRS) & ! LCOV_EXCL_LINE
                                                )
#endif
            if ( any(StateNew(1:nd)<=mc_DomainLowerLimitVec) .or. any(StateNew(1:nd)>=mc_DomainUpperLimitVec) ) then
                domainCheckCounter = domainCheckCounter + 1
                if (domainCheckCounter==mc_MaxNumDomainCheckToWarn) then
                    call warn( prefix = mc_methodBrand, outputUnit = self%logFileUnit, msg = mc_MaxNumDomainCheckToWarnMsg )
                end if
                if (domainCheckCounter==mc_MaxNumDomainCheckToStop) then
                    self%Err%occurred = .true.
                    self%Err%msg = mc_MaxNumDomainCheckToStopMsg
#if !defined CODECOV_ENABLED && ((!defined MATLAB_ENABLED && !defined PYTHON_ENABLED && !defined R_ENABLED) || defined CAF_ENABLED && defined MPI_ENABLED )
                    call abort( Err = self%Err, prefix = mc_methodBrand, newline = "\n", outputUnit = self%logFileUnit )
#endif
                    return
                end if
                cycle loopBoundaryCheck
            end if
            exit loopBoundaryCheck
        end do loopBoundaryCheck

    end subroutine getNew

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined PARADISE
    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> Return the natural logarithm of the probability of acceptance of the new state given the old state.
    !>
    !> @param[in]   nd          :   The number of dimensions of the objective function.
    !> @param[in]   counterDRS  :   The Delayed Rejection Stage counter.
    !> @param[in]   StateOld    :   The old sampled state.
    !> @param[in]   StateNew    :   The new sampled state.
    !>
    !> \return
    !> `logProb` : The log probability of obtaining obtaining the new sample given the old sample.
    ! LCOV_EXCL_START
    pure function getLogProb( nd                &
                            , counterDRS        &
                            , StateOld          &
                            , StateNew          &
                            ) result(logProb)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProb
#endif
#if defined NORMAL
        use Statistics_mod, only: getLogProbMVN
#elif defined UNIFORM
        use Statistics_mod, only: isInsideEllipsoid
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        implicit none
        integer(IK), intent(in)             :: nd
        integer(IK), intent(in)             :: counterDRS
        real(RK)   , intent(in)             :: StateOld(nd)
        real(RK)   , intent(in)             :: StateNew(nd)
        real(RK)                            :: logProb
#if defined UNIFORM
            if (isInsideEllipsoid(nd,StateNew-StateOld,mv_InvCovMat(1:mc_ndim,1:mc_ndim,counterDRS))) then
                logProb = mc_negLogVolUnitBall + mv_logSqrtDetInvCovMat(counterDRS)
            else
                logProb = NEGINF_RK
            end if
#elif defined NORMAL
            logProb = getLogProbMVN ( nd = nd &
                                    , MeanVec = StateOld &
                                    , InvCovMat = mv_InvCovMat(1:nd,1:nd,counterDRS) &
                                    , logSqrtDetInvCovMat = mv_logSqrtDetInvCovMat(counterDRS) &
                                    , Point = StateNew &
                                    )
#endif
    end function getLogProb
    ! LCOV_EXCL_STOP
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> Perform adaptation of the proposal distribution of the MCMC sampler.
    !>
    !> @param[in]       nd                      :   The number of dimensions of the objective function.
    !> @param[in]       chainSize               :   The size of the input chain in compact format.
    !> @param[in]       Chain                   :   The two dimensional `(nd, chainSize)` array of accepted states.
    !> @param[in]       ChainWeight             :   The one dimensional integer array of the weights of the sampled states in the chain.
    !> @param[in]       isFreshRun              :   The logical flag indicating whether the simulation is in regular (non-restart) mode.
    !> @param[in]       samplerUpdateIsGreedy   :   The logical flag indicating whether the adaptation is in greedy mode..
    !> @param[inout]    meanAccRateSinceStart   :   The mean acceptance rate since the start of the sampling.
    !> @param[out]      samplerUpdateSucceeded  :   The output logical flag indicating whether the sampling succeeded.
    !> @param[out]      adaptationMeasure       :   The output real number in the range `[0,1]` indicating the amount of adaptation,
    !>                                              with zero indicating no adaptation and one indicating extreme adaptation to the extent
    !>                                              that the new adapted proposal distribution is completely different from the previous proposal.
    !> \warning
    !> This routine must be exclusively called by the leader images.
    !>
    !> \remark
    !> For information on the meaning of `adaptationMeasure`, see the paper by Shahmoradi and Bagheri, 2020, whose PDF is available at:
    !> [https://www.cdslab.org/paramonte/notes/overview/preface/#the-paradram-sampler](https://www.cdslab.org/paramonte/notes/overview/preface/#the-paradram-sampler)
    subroutine doAdaptation ( self                      &
                            , nd                        &
                            , chainSize                 &
                            , Chain                     &
                            , ChainWeight               &
                            , isFreshRun                &
                            , samplerUpdateIsGreedy     &
                            , meanAccRateSinceStart     &
                            , samplerUpdateSucceeded    &
                            , adaptationMeasure         &
                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: doAdaptation
#endif

        use Statistics_mod, only: getSamCovUpperMeanTrans, getWeiSamCovUppMeanTrans, mergeMeanCovUpper ! LCOV_EXCL_LINE
        use Matrix_mod, only: getCholeskyFactor, getLogSqrtDetPosDefMat
        use Constants_mod, only: RK, IK ! , EPS_RK
        use String_mod, only: num2str
        use Err_mod, only: Err_type, abort, warn
        implicit none

        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME // "@doAdaptation()"

        class(Proposal_type), intent(inout)             :: self
        integer(IK) , intent(in)                        :: nd
        integer(IK) , intent(in)                        :: chainSize
        real(RK)    , intent(in)                        :: Chain(nd,chainSize)
        integer(IK) , intent(in)                        :: ChainWeight(chainSize)
        logical     , intent(in)                        :: isFreshRun
        logical     , intent(in)                        :: samplerUpdateIsGreedy
        real(RK)    , intent(inout)                     :: meanAccRateSinceStart ! is intent(out) in restart mode, intent(in) in fresh mode.
        logical     , intent(out)                       :: samplerUpdateSucceeded
        real(RK)    , intent(out)                       :: adaptationMeasure

        integer(IK)                                     :: i, j
        real(RK)                                        :: MeanNew(nd)
        real(RK)                                        :: MeanCurrent(nd)
        real(RK)                                        :: CovMatUpperOld(nd,nd)
        real(RK)                                        :: CovMatUpperCurrent(nd,nd)
        real(RK)                                        :: logSqrtDetNew
        real(RK)                                        :: logSqrtDetSum
        real(RK)                                        :: adaptiveScaleFactor
        logical                                         :: adaptationMeasureComputationNeeded ! only used to avoid redundant affinity computation, if no update occurs
        logical                                         :: singularityOccurred, scalingNeeded
        integer(IK)                                     :: sampleSizeOld, sampleSizeCurrent

        scalingNeeded = .false.
        sampleSizeOld = self%sampleSizeOld ! this is kept only for restoration of self%sampleSizeOld, if needed.

        ! read/write meanAccRateSinceStart from/to restart file

        if (mc_Image%isLeader) then
            if (.not. isFreshRun) then
                if (mc_isBinaryRestartFileFormat) then
                    call readRestartFileBinary(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = meanAccRateSinceStart)
                else
                    call readRestartFileAscii(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = meanAccRateSinceStart)
                end if
            end if
        end if

        ! First if there are less than nd+1 points for new covariance computation, then just scale the covariance and return

        blockSufficientSampleSizeCheck: if (chainSize>nd) then

            ! get the upper covariance matrix and Mean of the new sample

            if (samplerUpdateIsGreedy) then
                sampleSizeCurrent = chainSize
                call getSamCovUpperMeanTrans(sampleSizeCurrent,nd,Chain,CovMatUpperCurrent,MeanCurrent)
            else
                sampleSizeCurrent = sum(ChainWeight)
                call getWeiSamCovUppMeanTrans(chainSize,sampleSizeCurrent,nd,Chain,ChainWeight,CovMatUpperCurrent,MeanCurrent)
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! combine old and new covariance matrices if both exist

            blockMergeCovMat: if (self%sampleSizeOld==1_IK) then

                ! There is no prior old Covariance matrix to combine with the new one from the new chain

                self%MeanOld(1:nd) = MeanCurrent
                self%sampleSizeOld = sampleSizeCurrent

                ! copy and then scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                do j = 1, nd
                    do i = 1, j
                        CovMatUpperOld(i,j) = comv_CholDiagLower(i,j,0)   ! This will be used to recover the old covariance in case of update failure, and to compute the adaptation measure
                        comv_CholDiagLower(i,j,0) = CovMatUpperCurrent(i,j) * mc_defaultScaleFactorSq
                    end do
                end do

            else blockMergeCovMat

                ! first scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                do j = 1, nd
                    do i = 1, j
                        CovMatUpperOld(i,j) = comv_CholDiagLower(i,j,0)   ! This will be used to recover the old covariance in case of update failure, and to compute the adaptation measure
                        CovMatUpperCurrent(i,j) = CovMatUpperCurrent(i,j) * mc_defaultScaleFactorSq
                    end do
                end do

                ! now combine it with the old covariance matrix.
                ! Do not set the full boundaries' range `(1:nd)` for `comv_CholDiagLower` in the following subroutine call.
                ! Setting the boundaries forces the compiler to generate a temporary array.

                call mergeMeanCovUpper  ( nd            = nd                            & ! LCOV_EXCL_LINE
                                        , npA           = self%sampleSizeOld            & ! LCOV_EXCL_LINE
                                        , MeanVecA      = self%MeanOld                  & ! LCOV_EXCL_LINE
                                        , CovMatUpperA  = CovMatUpperOld                & ! LCOV_EXCL_LINE
                                        , npB           = sampleSizeCurrent             & ! LCOV_EXCL_LINE
                                        , MeanVecB      = MeanCurrent                   & ! LCOV_EXCL_LINE
                                        , CovMatUpperB  = CovMatUpperCurrent            & ! LCOV_EXCL_LINE
                                        , MeanVecAB     = MeanNew                       & ! LCOV_EXCL_LINE
                                        , CovMatUpperAB = comv_CholDiagLower(:,1:nd,0)  & ! LCOV_EXCL_LINE
                                        )
                self%MeanOld(1:nd) = MeanNew

            end if blockMergeCovMat

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! now get the Cholesky factorization

            ! WARNING: Do not set the full boundaries' range `(1:nd)` for the first index of `comv_CholDiagLower` in the following subroutine call.
            ! WARNING: Setting the boundaries forces the compiler to generate a temporary array.

            call getCholeskyFactor( nd, comv_CholDiagLower(:,1:nd,0), comv_CholDiagLower(1:nd,0,0) )

            blockPosDefCheck: if (comv_CholDiagLower(1,0,0)>0._RK) then

                !singularityOccurred = .false.
                samplerUpdateSucceeded = .true.
                adaptationMeasureComputationNeeded = .true.
                self%sampleSizeOld = self%sampleSizeOld + sampleSizeCurrent

            ! LCOV_EXCL_START
            else blockPosDefCheck

                adaptationMeasure = 0._RK
                !singularityOccurred = .true.
                samplerUpdateSucceeded = .false.
                adaptationMeasureComputationNeeded = .false.
                self%sampleSizeOld = sampleSizeOld

                ! it may be a good idea to add a warning message printed out here for the singularity occurrence

                call warn   ( prefix = mc_methodBrand       & ! LCOV_EXCL_LINE
                            , outputUnit = self%logFileUnit & ! LCOV_EXCL_LINE
                            , marginTop = 0_IK              & ! LCOV_EXCL_LINE
                            , marginBot = 0_IK              & ! LCOV_EXCL_LINE
                            , msg = "Singularity occurred while updating the proposal distribution's covariance matrix." & ! LCOV_EXCL_LINE
                            )

                ! recover the old upper covariance matrix

                do j = 1, nd
                    do i = 1, j
                        comv_CholDiagLower(i,j,0) = CovMatUpperOld(i,j)
                    end do
                end do

                ! ensure the old Cholesky factorization can be recovered

                ! WARNING: Do not set the full boundaries' range `(1:nd)` for the first index of `comv_CholDiagLower` in the following subroutine call.
                ! WARNING: Setting the boundaries forces the compiler to generate a temporary array.

                call getCholeskyFactor( nd, comv_CholDiagLower(:,1:nd,0), comv_CholDiagLower(1:nd,0,0) ) ! avoid temporary array creation by using : indexer
                if (comv_CholDiagLower(1,0,0)<0._RK) then
                    write(self%logFileUnit,"(A)")
                    write(self%logFileUnit,"(A)") "Singular covariance matrix detected:"
                    write(self%logFileUnit,"(A)")
                    do j = 1, nd
                        write(self%logFileUnit,"(*(E25.15))") comv_CholDiagLower(1:j,j,0)
                    end do
                    write(self%logFileUnit,"(A)")
                    self%Err%occurred = .true.
                    self%Err%msg =  PROCEDURE_NAME // &
                                    ": Error occurred while attempting to compute the Cholesky factorization of the &
                                    &covariance matrix of the proposal distribution of " // mc_methodName // "'s sampler. &
                                    &This is highly unusual, and can be indicative of some major underlying problems.\n&
                                    &It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations. &
                                    &In such case, consider increasing the value of the input variable adaptiveUpdatePeriod.\n&
                                    &It may also be that your input objective function has been incorrectly implemented.\n&
                                    &For example, ensure that you are passing a correct value of ndim to the ParaMonte sampler routine,\n&
                                    &the same value that is expected as input to your objective function's implementation.\n&
                                    &Otherwise, restarting the simulation might resolve the error."
                    call abort( Err = self%Err, prefix = mc_methodBrand, newline = "\n", outputUnit = self%logFileUnit )
                    return
                end if

            end if blockPosDefCheck
            ! LCOV_EXCL_STOP

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !! perform global adaptive scaling is requested
            !if (mc_scalingRequested) then
            !    if (meanAccRateSinceStart<mc_targetAcceptanceRate) then
            !        self%adaptiveScaleFactorSq = mc_maxScaleFactorSq**((mc_targetAcceptanceRate-meanAccRateSinceStart)/mc_targetAcceptanceRate)
            !    else
            !        self%adaptiveScaleFactorSq = mc_maxScaleFactorSq**((mc_targetAcceptanceRate-meanAccRateSinceStart)/(1._RK-mc_targetAcceptanceRate))
            !    end if
            !    self%adaptiveScaleFactorSq = self%adaptiveScaleFactorSq * (meanAccRateSinceStart/mc_targetAcceptanceRate)**mc_ndimInverse
            !end if
            if (mc_scalingRequested) scalingNeeded = .true.

        else blockSufficientSampleSizeCheck

            ! singularity has occurred. If the first covariance merging has not occurred yet, set the scaling factor appropriately to shrink the covariance matrix.

            samplerUpdateSucceeded = .false.
            if (self%sampleSizeOld==1_IK .or. mc_scalingRequested) then
                scalingNeeded = .true.
                adaptationMeasureComputationNeeded = .true.
                ! save the old covariance matrix for the computation of the adaptation measure
                do j = 1, nd
                    do i = 1,j
                        CovMatUpperOld(i,j) = comv_CholDiagLower(i,j,0)
                    end do
                end do
            else
                adaptationMeasure = 0._RK
                adaptationMeasureComputationNeeded = .false.
            end if


        end if blockSufficientSampleSizeCheck

        ! adjust the scale of the covariance matrix and the Cholesky factor, if needed

!scalingNeeded = .true.
        if (scalingNeeded) then
            if ( meanAccRateSinceStart < mc_TargetAcceptanceRateLimit(1) .or. meanAccRateSinceStart > mc_TargetAcceptanceRateLimit(2) ) then
                self%adaptiveScaleFactorSq = (meanAccRateSinceStart/mc_targetAcceptanceRate)**mc_ndimInverse
!block
!    use Statistics_mod, only: getRandUniform
!    integer, save :: counter = 0_IK
!    counter = counter - 1
!    self%adaptiveScaleFactorSq = self%adaptiveScaleFactorSq * exp(-counter*getRandUniform(-1.e0_RK,1.e0_RK)/1.e4_RK)
!    !use Statistics_mod, only: getRandInt
!    !self%adaptiveScaleFactorSq = self%adaptiveScaleFactorSq * exp(real(getRandInt(-1_IK,1_IK),kind=RK))
!    !write(*,*) counter, self%adaptiveScaleFactorSq
!end block
                adaptiveScaleFactor = sqrt(self%adaptiveScaleFactorSq)
                do j = 1, nd
                    ! update the Cholesky diagonal elements
                    comv_CholDiagLower(j,0,0) = comv_CholDiagLower(j,0,0) * adaptiveScaleFactor
                    ! update covariance matrix
                    do i = 1,j
                        comv_CholDiagLower(i,j,0) = comv_CholDiagLower(i,j,0) * self%adaptiveScaleFactorSq
                    end do
                    ! update the Cholesky factorization
                    do i = j+1, nd
                        comv_CholDiagLower(i,j,0) = comv_CholDiagLower(i,j,0) * adaptiveScaleFactor
                    end do
                end do
            end if
        end if

        ! compute the adaptivity only if any updates has occurred


        blockAdaptationMeasureComputation: if (adaptationMeasureComputationNeeded) then

            logSqrtDetNew = sum(log( comv_CholDiagLower(1:nd,0,0) ))

            ! use the universal upper bound

            do j = 1, nd
                do i = 1,j
                    CovMatUpperCurrent(i,j) = 0.5_RK * ( comv_CholDiagLower(i,j,0) + CovMatUpperOld(i,j) ) ! dummy
                end do
            end do
            call getLogSqrtDetPosDefMat(nd,CovMatUpperCurrent,logSqrtDetSum,singularityOccurred)
            if (singularityOccurred) then
                ! LCOV_EXCL_START
                write(self%logFileUnit,"(A)")
                write(self%logFileUnit,"(A)") "Singular covariance matrix detected while computing the Adaptation measure:"
                write(self%logFileUnit,"(A)")
                do j = 1, nd
                    write(self%logFileUnit,"(*(E25.15))") CovMatUpperCurrent(1:j,j)
                end do
                self%Err%occurred = .true.
                self%Err%msg =  PROCEDURE_NAME // &
                                ": Error occurred while computing the Cholesky factorization of &
                                &a matrix needed for the computation of the Adaptation measure. &
                                &Such error is highly unusual, and requires an in depth investigation of the case.\n&
                                &It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations. &
                                &In such case, consider increasing the value of the input variable adaptiveUpdatePeriod.\n&
                                &It may also be that your input objective function has been incorrectly implemented.\n&
                                &For example, ensure that you are passing a correct value of ndim to the ParaMonte sampler routine,\n&
                                &the same value that is expected as input to your objective function's implementation.\n&
                                &Otherwise, restarting the simulation might resolve the error."
                call abort( Err = self%Err, prefix = mc_methodBrand, newline = "\n", outputUnit = self%logFileUnit )
                return
                ! LCOV_EXCL_STOP
            end if

            !adaptationMeasure = 1._RK - exp( 0.5_RK*(self%logSqrtDetOld+logSqrtDetNew) - logSqrtDetSum )
            adaptationMeasure = 1._RK - exp( 0.5*(self%logSqrtDetOld + logSqrtDetNew) - logSqrtDetSum ) ! totalVariationUpperBound
            if (adaptationMeasure>0._RK) then
                adaptationMeasure = sqrt(adaptationMeasure) ! totalVariationUpperBound
            ! LCOV_EXCL_START
            elseif (adaptationMeasure<0._RK) then
                call warn   ( prefix = mc_methodBrand &
                            , outputUnit = self%logFileUnit &
                            , msg = mc_negativeTotalVariationMsg//num2str(adaptationMeasure) )
                adaptationMeasure = 0._RK
            end if
            ! LCOV_EXCL_STOP
            self%logSqrtDetOld = logSqrtDetNew

!block
!integer, save :: counter = 0
!counter = counter + 1
!!if (counter==1) then
!if (adaptationMeasure>1._RK) then
!write(*,*)
!write(*,*) self%logSqrtDetOld
!write(*,*) logSqrtDetNew
!write(*,*) logSqrtDetSum
!write(*,*) self%logSqrtDetOld + logSqrtDetNew - 2_IK * logSqrtDetSum
!write(*,*) exp( self%logSqrtDetOld + logSqrtDetNew - 2_IK * logSqrtDetSum )
!write(*,*)
!end if
!end block

            ! update the higher-stage delayed-rejection Cholesky Lower matrices

            if (mc_delayedRejectionRequested) call updateDelRejCholDiagLower()
            call getInvCovMat()

        end if blockAdaptationMeasureComputation

        ! read/write restart file

        if (mc_Image%isLeader) then
            if (isFreshRun) then
                if (mc_isBinaryRestartFileFormat) then
                    call writeRestartFileBinary(restartFileUnit = self%restartFileUnit, meanAccRateSinceStart = meanAccRateSinceStart)
                elseif (mc_isAsciiRestartFileFormat) then
                    call writeRestartFileAscii  ( restartFileUnit = self%restartFileUnit & ! LCOV_EXCL_LINE
                                                , meanAccRateSinceStart = meanAccRateSinceStart & ! LCOV_EXCL_LINE
                                                , sampleSizeOld = self%sampleSizeOld & ! LCOV_EXCL_LINE
                                                , logSqrtDetOld = self%logSqrtDetOld & ! LCOV_EXCL_LINE
                                                , adaptiveScaleFactorSq = self%adaptiveScaleFactorSq & ! LCOV_EXCL_LINE
                                                , MeanOld = self%MeanOld & ! LCOV_EXCL_LINE
                                                )
                end if
            end if
        end if

    end subroutine doAdaptation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! LCOV_EXCL_START
    ! ATTN: This routine needs further correction for the delayed rejection method
#if false
    subroutine doAutoTune   ( adaptationMeasure &
                            , AutoTuneScaleSq   &
                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: doAutoTune
#endif
        use Matrix_mod, only: getLogSqrtDetPosDefMat
        use Constants_mod, only: RK, IK
        use Err_mod, only: abort
        implicit none
        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME // "@doAutoTune()"

        real(RK)   , intent(in)     :: AutoTuneScaleSq(1)
        real(RK)   , intent(inout)  :: adaptationMeasure
        real(RK)                    :: logSqrtDetSum, logSqrtDetOld, logSqrtDetNew
        real(RK)                    :: CovMatUpperOld(1,1), CovMatUpperCurrent(1,1)
        logical                     :: singularityOccurred

        CovMatUpperOld = comv_CholDiagLower(1:mc_ndim,1:mc_ndim,0)
        logSqrtDetOld = sum(log( comv_CholDiagLower(1:mc_ndim,0,0) ))

        if (AutoTuneScaleSq(1)==0._RK) then
            comv_CholDiagLower(1,1,0) = 0.25_RK*comv_CholDiagLower(1,1,0)
            comv_CholDiagLower(1,0,0) = sqrt(comv_CholDiagLower(1,1,0))
        else
            comv_CholDiagLower(1,1,0) = AutoTuneScaleSq(1)
            comv_CholDiagLower(1,0,0) = sqrt(AutoTuneScaleSq(1))
        end if

        ! compute the adaptivity

        logSqrtDetNew = sum(log( comv_CholDiagLower(1:mc_ndim,0,0) ))
        CovMatUpperCurrent = 0.5_RK * ( comv_CholDiagLower(1:mc_ndim,1:mc_ndim,0) + CovMatUpperOld )
        call getLogSqrtDetPosDefMat(1_IK,CovMatUpperCurrent,logSqrtDetSum,singularityOccurred)
        if (singularityOccurred) then
            self%Err%occurred = .true.
            self%Err%msg =  PROCEDURE_NAME // &
                            ": Error occurred while computing the Cholesky factorization of &
                            &a matrix needed for the computation of the proposal distribution's adaptation measure. &
                            &Such error is highly unusual, and requires an in depth investigation of the case. &
                            &It may also be that your input objective function has been incorrectly implemented.\n&
                            &For example, ensure that you are passing a correct value of ndim to the ParaMonte sampler routine,\n&
                            &the same value that is expected as input to your objective function's implementation.\n&
                            &Otherwise, restarting the simulation might resolve the error."
            call abort( Err = self%Err, prefix = mc_methodBrand, newline = "\n", outputUnit = self%logFileUnit )
            return
        end if
        adaptationMeasure = 1._RK - exp( 0.5_RK*(logSqrtDetOld+logSqrtDetNew) - logSqrtDetSum )

    end subroutine doAutoTune
    ! LCOV_EXCL_STOP
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Note: based on some benchmarks with ndim = 1, the new design with merging cholesky diag and lower is faster than the original
    ! Note: double-communication, implementation. Here are some timings on 4 images:
    ! Note: new single communication:
    ! Note: image 2: avgTime =  6.960734060198531E-006
    ! Note: image 3: avgTime =  7.658279491640721E-006
    ! Note: image 4: avgTime =  9.261191328273417E-006
    ! Note: avg(avgTime): 7.960068293370891e-06
    ! Note: old double communication:
    ! Note: image 4: avgTime =  1.733505615104837E-005
    ! Note: image 3: avgTime =  1.442608268140151E-005
    ! Note: image 2: avgTime =  1.420345299036357E-005
    ! Note: avg(avgTime): 1.532153060760448e-05
    ! Note: avg(speedup): 1.924798889020109
    ! Note: One would expect this speed up to diminish as ndim goes to infinity,
    ! Note: since data transfer will dominate the communication overhead.
    !> \brief
    !> Broadcast adaptation to all images.
    !> \warning
    !> When CAF parallelism is used, this routine must be exclusively called by the rooter images.
    !> When MPI parallelism is used, this routine must be called by all images.
#if defined CAF_ENABLED || MPI_ENABLED
    subroutine bcastAdaptation()
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: bcastAdaptation
#endif
#if defined CAF_ENABLED
        implicit none
        comv_CholDiagLower(1:mc_ndim,0:mc_ndim,0) = comv_CholDiagLower(1:mc_ndim,0:mc_ndim,0)[1]
        if (mc_delayedRejectionRequested) call updateDelRejCholDiagLower()  ! update the higher-stage delayed-rejection Cholesky Lower matrices
#elif defined MPI_ENABLED
        use mpi ! LCOV_EXCL_LINE
        implicit none
        integer :: ierrMPI
        call mpi_bcast  ( comv_CholDiagLower    & ! LCOV_EXCL_LINE ! buffer: XXX: first element need not be shared. This may need a fix in future.
                        , mc_ndimSqPlusNdim     & ! LCOV_EXCL_LINE ! count
                        , mpi_double_precision  & ! LCOV_EXCL_LINE ! datatype
                        , 0                     & ! LCOV_EXCL_LINE ! root: broadcasting rank
                        , mpi_comm_world        & ! LCOV_EXCL_LINE ! comm
                        , ierrMPI               & ! LCOV_EXCL_LINE ! ierr
                        )
        ! It is essential for the following to be exclusively done by the rooter images. The leaders have had their updates in `doAdaptation()`.
        if (mc_Image%isRooter .and. mc_delayedRejectionRequested) call updateDelRejCholDiagLower()
#endif
        call getInvCovMat()
    end subroutine bcastAdaptation
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse covariance matrix of the current covariance of the proposal distribution.
    !>
    !> \warning
    !> This procedure is `impure` and sets the values of multiple global variables upon exit.
    !>
    !> \remark
    !> This procedure is required for the proposal probabilities
    !>
    !> \todo
    !> The specification of the array bounds causes the compiler to generate a temporary copy of the array, despite being unnecessary.
    !> Right now, this is resolved by replacing the array bounds with `:`. A better solution is to add the `contiguous` attribute to
    !> the corresponding argument of [getInvMatFromCholFac](ref matrix_mod::getinvmatfromcholfac) to guarantee it to the compiler.
    !> More than improving performance, this would turn off the pesky compiler warnings about temporary array creation.
    subroutine getInvCovMat()
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvCovMat
#endif
        use Matrix_mod, only: getInvMatFromCholFac ! LCOV_EXCL_LINE
        implicit none
        integer(IK) :: istage
        ! update the inverse covariance matrix of the proposal from the computed Cholesky factor
        do concurrent(istage=0:mc_DelayedRejectionCount)
            ! WARNING: Do not set the full boundaries' range `(1:mc_ndim)` for the first index of `comv_CholDiagLower` in the following subroutine call.
            ! WARNING: Setting the boundaries forces the compiler to generate a temporary array.
            mv_InvCovMat(1:mc_ndim,1:mc_ndim,istage) = getInvMatFromCholFac ( nd = mc_ndim & ! LCOV_EXCL_LINE
                                                                            , CholeskyLower = comv_CholDiagLower(:,1:mc_ndim,istage) & ! LCOV_EXCL_LINE
                                                                            , CholeskyDiago = comv_CholDiagLower(1:mc_ndim,0,istage) & ! LCOV_EXCL_LINE
                                                                            )
            mv_logSqrtDetInvCovMat(istage) = -sum(log( comv_CholDiagLower(1:mc_ndim,0,istage) ))
        end do
!if (mc_ndim==1 .and. abs(log(sqrt(mv_InvCovMat(1,1,0)))-mv_logSqrtDetInvCovMat(0))>1.e-13_RK) then
!write(*,"(*(g0,:,' '))") "log(sqrt(mv_InvCovMat(1,1,0))) /= mv_logSqrtDetInvCovMat(0)"
!write(*,"(*(g0,:,' '))") log(sqrt(mv_InvCovMat(1,1,0))), mv_logSqrtDetInvCovMat(0)
!write(*,"(*(g0,:,' '))") abs(-log(sqrt(mv_InvCovMat(1,1,0)))-mv_logSqrtDetInvCovMat(0))
!error stop
!endif
    end subroutine getInvCovMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Update the Cholesky factorization of the higher-stage delayed rejection covariance matrices.
    !>
    !> \warning
    !> This procedure is impure and sets the values of multiple global variables upon exit.
    !>
    !> \todo
    !> The performance of this update could be improved by only updating the higher-stage covariance, only when needed.
    !> However, the gain will be likely minimal, especially in low-dimensions.
    subroutine updateDelRejCholDiagLower()
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: updateDelRejCholDiagLower
#endif
        implicit none
        integer(IK) :: j, istage
        ! update the Cholesky factor of the delayed-rejection-stage proposal distributions
        do istage = 1, mc_DelayedRejectionCount
            comv_CholDiagLower(1:mc_ndim,0,istage) = comv_CholDiagLower(1:mc_ndim,0,istage-1) * mc_DelayedRejectionScaleFactorVec(istage)
            do j = 1, mc_ndim
                comv_CholDiagLower(j+1:mc_ndim,j,istage) = comv_CholDiagLower(j+1:mc_ndim,j,istage-1) * mc_DelayedRejectionScaleFactorVec(istage)
            end do
        end do
        ! There is no need to check for positive-definiteness of the comv_CholDiagLower, it is already checked on the first image.
    end subroutine updateDelRejCholDiagLower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> This procedure is called by the sampler kernel routines.\n
    !> Write the restart information to the output file.
    !>
    !> @param[in]   meanAccRateSinceStart   : The current mean acceptance rate of the sampling (**optional**).
    !> @param[in]   MeanOld                 : The mean of the old sample (**optional**, to be present only if `meanAccRateSinceStart` is missing.).
    !>
    !> \warning
    !> The input argument `MeanOld` must be present if and only if `meanAccRateSinceStart` is missing as an input arguments.
    !> This condition will **NOT** be checked for at runtime. It is the developer's responsibility to ensure it holds.
    subroutine writeRestartFileAscii(restartFileUnit, meanAccRateSinceStart, sampleSizeOld, logSqrtDetOld, adaptiveScaleFactorSq, MeanOld)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writeRestartFileAscii
#endif
        implicit none
        integer(IK) , intent(in)    :: restartFileUnit
        real(RK)    , intent(in)    :: meanAccRateSinceStart
        integer(IK) , intent(in)    :: sampleSizeOld
        real(RK)    , intent(in)    :: logSqrtDetOld
        real(RK)    , intent(in)    :: adaptiveScaleFactorSq
        real(RK)    , intent(in)    :: MeanOld(mc_ndim)
        integer(IK)                 :: i, j
        write(restartFileUnit, mc_restartFileFormat ) "meanAcceptanceRateSinceStart", meanAccRateSinceStart
        write(restartFileUnit, mc_restartFileFormat ) "sampleSize" &
                                                    , sampleSizeOld &
                                                    , "logSqrtDeterminant" &
                                                    , logSqrtDetOld &
                                                    , "adaptiveScaleFactorSquared" &
                                                    , adaptiveScaleFactorSq * mc_defaultScaleFactorSq &
                                                    , "meanVec" &
                                                    , MeanOld(1:mc_ndim) & ! LCOV_EXCL_LINE
                                                    , "covMat" &
                                                    , ((comv_CholDiagLower(i,j,0),i=1,j),j=1,mc_ndim)
                                                    !, (comv_CholDiagLower(1:mc_ndim,0:mc_ndim,0)
        flush(restartFileUnit)
    end subroutine writeRestartFileAscii

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> This procedure is called by the sampler kernel routines.\n
    !> Write the restart information to the output file.
    !>
    !> @param[in]   meanAccRateSinceStart   : The current mean acceptance rate of the sampling (**optional**).
    !> @param[in]   MeanOld                 : The mean of the old sample (**optional**, to be present only if `meanAccRateSinceStart` is missing.).
    !>
    !> \warning
    !> The input argument `MeanOld` must be present if and only if `meanAccRateSinceStart` is missing as an input arguments.
    !> This condition will **NOT** be checked for at runtime. It is the developer's responsibility to ensure it holds.
    subroutine writeRestartFileBinary(restartFileUnit, meanAccRateSinceStart)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writeRestartFileBinary
#endif
        implicit none
        integer(IK) , intent(in) :: restartFileUnit
        real(RK)    , intent(in) :: meanAccRateSinceStart
        write(restartFileUnit) meanAccRateSinceStart
        flush(restartFileUnit)
    end subroutine writeRestartFileBinary

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> This procedure is called by the sampler kernel routines.\n
    !> Read the restart information from the restart file.
    !>
    !> @param[out]  meanAccRateSinceStart : The current mean acceptance rate of the sampling (**optional**).
    subroutine readRestartFileAscii(restartFileUnit, meanAccRateSinceStart)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: readRestartFileAscii
#endif
        implicit none
        integer(IK) , intent(in)    :: restartFileUnit
        real(RK)    , intent(out)   :: meanAccRateSinceStart
        integer(IK)                 :: i
        read(restartFileUnit,*)
        read(restartFileUnit,*) meanAccRateSinceStart
        do i = 1, 8 + mc_ndim * (mc_ndim+3) / 2
            read(restartFileUnit, *)
        end do
    end subroutine readRestartFileAscii

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the [ParaXXXX_ProposalNormal_type](@ref paraxxxxproposalnormal_type)
    !> or [ParaXXXX_ProposalUniform_type](@ref paraxxxxproposaluniform_type) classes.\n
    !> This procedure is called by the sampler kernel routines.\n
    !> Read the restart information from the restart file.
    !>
    !> @param[out]  meanAccRateSinceStart : The current mean acceptance rate of the sampling (**optional**).
    subroutine readRestartFileBinary(restartFileUnit, meanAccRateSinceStart)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: readRestartFileBinary
#endif
        implicit none
        integer(IK) , intent(in)    :: restartFileUnit
        real(RK)    , intent(out)   :: meanAccRateSinceStart
        read(restartFileUnit) meanAccRateSinceStart
    end subroutine readRestartFileBinary

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#undef ParaXXXX
#undef GET_RANDOM_PROPOSAL

