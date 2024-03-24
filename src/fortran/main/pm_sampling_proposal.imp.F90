!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if     CAF_ENABLED
#define SET_CAF(X)X
#else
#define SET_CAF(X)
#endif
        ! \bug
        ! Avoid Intel ifort bug for too many `use` statements in a submodule by placing them all in the submodule header.
#if     CHECK_ENABLED
        use pm_err, only: setAsserted, getFine
#define CHECK_ASSERTION(LINE,ASSERTION,MSG)call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
        ! Define error handling.
#define RETURN_IF_FAILED(LINE,FAILED,MSG) \
if (FAILED) then; \
err%occurred = .true._LK; \
err%msg = PROCEDURE_NAME//getLine(LINE)//SKC_": "//MSG; \
call spec%disp%stop%show(err%msg); \
return; \
end if;
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     ParaDISE_ENABLED || ParaDRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_err, only: getLine
        use pm_err, only: err_type
        use pm_val2str, only: getStr
        use pm_kind, only: SKC => SK, SK, IK, LK
        use pm_paramonte, only: PARAMONTE_WEB_ISSUES
#if     ParaDISE_ENABLED
        use pm_sampling_dise, only: NL1, NL2
        use pm_sampling_dise, only: spec_type => specdise_type
        use pm_sampling_dise, only: stat_type => statdise_type
        use pm_sampling_scio, only: getErrChainRead => getErrChainReadDISE
        use pm_sampling_scio, only: getErrChainWrite => getErrChainWriteDISE
        use pm_sampling_scio, only: chainFileColName => chainFileColNameDISE
        use pm_sampling_scio, only: isFailedChainResize => isFailedChainResizeDISE
        implicit none
        character(*,SKC), parameter :: MODULE_NAME = SK_"@pm_sampling_proposal_dise"
#elif   ParaDRAM_ENABLED
        use pm_sampling_dram, only: NL1, NL2
        use pm_sampling_dram, only: spec_type => specdram_type
        use pm_sampling_dram, only: stat_type => statdram_type
        use pm_sampling_scio, only: getErrChainRead => getErrChainReadDRAM
        use pm_sampling_scio, only: getErrChainWrite => getErrChainWriteDRAM
        use pm_sampling_scio, only: chainFileColName => chainFileColNameDRAM
        use pm_sampling_scio, only: isFailedChainResize => isFailedChainResizeDRAM
        implicit none
        character(*,SKC), parameter :: MODULE_NAME = SK_"@pm_sampling_proposal_dram"
#endif
#if     CAF_ENABLED
        real(RKC), save     , allocatable   :: comv_covLowChoUpp(:,:)[:] ! lower covariance, upper Cholesky.
#endif
        integer(IK)         , private       :: ndimp1
        integer(IK)         , private       :: ndim
        type :: scaleFactorSq_type
            real(RKC) :: default
            real(RKC) :: running
        end type
        type :: proposal_type
            ! the following are made components for the sake of thread-safe save attribute and the restart file generation.
            character(:,SKC), allocatable   :: format
            type(scaleFactorSq_type)        :: scaleFactorSq
            logical(LK)                     :: delRejEnabled
            integer(IK)                     :: sampleSizeOld
            real(RKC)                       :: logSqrtDetOld
            real(RKC)       , allocatable   :: covLowChoUpp(:,:,:) ! The covariance Matrix of the proposal distribution. Last index belongs to delayed rejection.
            real(RKC)       , allocatable   :: meanOld(:)
#if         ParaDISE_ENABLED
            real(RKC)                       :: negLogVolUnitBall
            real(RKC)       , allocatable   :: logSqrtDetInvCov(:) ! (0 : spec%proposalDelayedRejectionCount%val)
            real(RKC)       , allocatable   :: invCov(:,:,:)
#define     SET_ParaDISE(X)X
#else
#define     SET_ParaDISE(X)
#endif
        end type
        interface proposal_type
            procedure :: constructProposal
        end interface

contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine killMeAlreadyCMake1_RK5(); use pm_sampling_scio_RK5, only: RKC; end subroutine
        subroutine killMeAlreadyCMake1_RK4(); use pm_sampling_scio_RK4, only: RKC; end subroutine
        subroutine killMeAlreadyCMake1_RK3(); use pm_sampling_scio_RK3, only: RKC; end subroutine
        subroutine killMeAlreadyCMake1_RK2(); use pm_sampling_scio_RK2, only: RKC; end subroutine
        subroutine killMeAlreadyCMake1_RK1(); use pm_sampling_scio_RK1, only: RKC; end subroutine
#if     ParaDISE_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_dise_RK5, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_dise_RK4, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_dise_RK3, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_dise_RK2, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_dise_RK1, only: RKC; end subroutine
#elif   ParaDRAM_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_dram_RK5, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_dram_RK4, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_dram_RK3, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_dram_RK2, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_dram_RK1, only: RKC; end subroutine
#elif   ParaNest_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_nest_RK5, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_nest_RK4, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_nest_RK3, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_nest_RK2, only: RKC; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_nest_RK1, only: RKC; end subroutine
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CAF_ENABLED || MPI_ENABLED
        !>  \brief
        !>  Broadcast adaptation to all images.
        !>
        !> \warning
        !> When CAF parallelism is used, this routine must be first called by the leader image, then exclusively called by the rooter images.
        !> When MPI parallelism is used, this routine must be called by all images.
        !>
        !>  \devnote
        !>  based on some benchmarks with ndim = 1, the new design with merging cholesky diag and lower is faster than the original
        !>  double-communication, implementation. Here are some timings on 4 images:
        !>  new single communication:
        !>  image 2: avgTime =  6.960734060198531E-006
        !>  image 3: avgTime =  7.658279491640721E-006
        !>  image 4: avgTime =  9.261191328273417E-006
        !>  avg(avgTime): 7.960068293370891e-06
        !>  old double communication:
        !>  image 4: avgTime =  1.733505615104837E-005
        !>  image 3: avgTime =  1.442608268140151E-005
        !>  image 2: avgTime =  1.420345299036357E-005
        !>  avg(avgTime): 1.532153060760448e-05
        !>  avg(speedup): 1.924798889020109
        !>  One would expect this speed up to diminish as ndim goes to infinity,
        !>  since data transfer will dominate the communication overhead.
        subroutine bcastProposalAdaptation(spec, proposal)
#if         CAF_ENABLED
            type(proposal_type), intent(inout) :: proposal
            type(spec_type), intent(in) :: spec
            if (spec%image%is%leader) then
                comv_covLowChoUpp(1 : ndim, 1 : ndimp1) = proposal%covLowChoUpp(1: ndim, 1 : ndimp1, 0)
            else
                proposal%covLowChoUpp(1: ndim, 1 : ndimp1, 0) = comv_covLowChoUpp(1 : ndim, 1 : ndimp1)[1]
                if (proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal) ! update the higher-stage delayed-rejection Cholesky Upper matrices.
                SET_ParaDISE(call setProposalInvCov(proposal))
            end if
#elif       MPI_ENABLED
            use mpi !mpi_f08, only: mpi_bcast, mpi_double_precision, mpi_comm_world ! LCOV_EXCL_LINE
            use, intrinsic :: iso_c_binding, only: c_sizeof
            integer, parameter :: NBRK = c_sizeof(0._RKC) ! Number of bytes in the real kind.
            type(spec_type), intent(in) :: spec
            type(proposal_type), intent(inout) :: proposal
            integer :: ierrMPI
            call mpi_bcast  ( proposal%covLowChoUpp & ! LCOV_EXCL_LINE
                           !, size(proposal%covLowChoUpp) & ! LCOV_EXCL_LINE ! count
                           !, mpi_double_precision & ! LCOV_EXCL_LINE ! datatype
                            , size(proposal%covLowChoUpp) * NBRK & ! LCOV_EXCL_LINE ! count
                            , mpi_byte & ! LCOV_EXCL_LINE ! datatype
                            , 0 & ! LCOV_EXCL_LINE ! root: broadcasting rank
                            , mpi_comm_world & ! LCOV_EXCL_LINE ! comm
                            , ierrMPI & ! LCOV_EXCL_LINE ! ierr
                            )
            ! It is essential for the following to be exclusively done by the rooter images. The leaders have had their updates in `setProposalAdapted()`.
            if (spec%image%is%rooter .and. proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal)
            SET_ParaDISE(call setProposalInvCov(spec, proposal))
#endif
        end subroutine bcastProposalAdaptation
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ParaDISE_ENABLED
        !>  \brief
        !>  Return the inverse covariance matrix of the current covariance of the proposal distribution.
        !>
        !>  \remark
        !>  This procedure is required for the proposal probabilities
        !>
        !>  \todo
        !>  The specification of the array bounds causes the compiler to generate a temporary copy of the array, despite being unnecessary.
        !>  Right now, this is resolved by replacing the array bounds with `:`. A better solution is to add the `contiguous` attribute to
        !>  the corresponding argument of [getMatInvFromChoLow](ref pm_matrixInv::getMatInvfromcholow) to guarantee it to the compiler.
        !>  More than improving performance, this would turn off the pesky compiler warnings about temporary array creation.
        PURE subroutine setProposalInvCov(spec, proposal)
            use pm_matrixInv, only: setMatInv, choUpp
            type(proposal_type), intent(inout) :: proposal
            type(spec_type), intent(in) :: spec
            integer(IK) :: istage
            ! update the inverse covariance matrix of the proposal from the computed Cholesky factor.
            istage = 0
            call setMatInv(proposal%invCov(:, 1 : ndim, istage), proposal%covLowChoUpp(:, 2 : ndimp1, istage), choUpp)
            proposal%logSqrtDetInvCov(istage) = -proposal%logSqrtDetOld
            do istage = 0, spec%proposalDelayedRejectionCount%val
                call setMatInv(proposal%invCov(:, 1 : ndim, istage), proposal%covLowChoUpp(:, 2 : ndimp1, istage), choUpp)
                !proposal%logSqrtDetInvCov(istage) = .5_RKC * getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : ndimp1, istage))
                proposal%logSqrtDetInvCov(istage) = proposal%logSqrtDetInvCov(istage) - ndim * log(spec%proposalDelayedRejectionScaleFactor%val(istage))
            end do
        end subroutine setProposalInvCov

        !>  \brief
        !>  Generate and return the natural logarithm of the probability of acceptance of the new state given the old state.
        PURE function getProposalJumpLogProb(spec, proposal, delayedRejectionStage, origin, destin) result(logPDF)
            use pm_distMultiNorm, only: getMultiNormLogPDFNF
            use pm_distMultiNorm, only: getMultiNormLogPDF
            use pm_ellipsoid, only: isMemberEll
            type(spec_type), intent(in) :: spec
            type(proposal_type), intent(in) :: proposal
            integer(IK), intent(in) :: delayedRejectionStage
            real(RKC), intent(in), contiguous :: origin(:), destin(:)
            character(*,SKC), parameter :: MSG = SKC_"@getProposalJumpLogProb()/: Internal library error; The specified uniform proposal jump is not within the proposal domain. "
            character(*,SKC), parameter :: ERRMSG = MSG//SKC_"ndim, delayedRejectionStage, proposal%invCov(:,1 : ndim, delayedRejectionStage), origin, destin = "
            real(RKC) :: logPDF
            if (spec%proposal%is%normal) then
                logPDF = getMultiNormLogPDF(destin, origin, proposal%invCov(:,1 : ndim, delayedRejectionStage), logPDFNF = getMultiNormLogPDFNF(ndim, proposal%logSqrtDetInvCov(delayedRejectionStage)))
            elseif (spec%proposal%is%uniform) then
                logPDF = proposal%negLogVolUnitBall + proposal%logSqrtDetInvCov(delayedRejectionStage)
                CHECK_ASSERTION(__LINE__, isMemberEll(proposal%invCov(:,1 : ndim, delayedRejectionStage), origin, destin), \
                ERRMSG//getStr([ndim, delayedRejectionStage])//SKC_", "//getStr([proposal%invCov(:,1 : ndim, delayedRejectionStage), origin, destin]))
            end if
        end function getProposalJumpLogProb
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function constructProposal(spec, err) result(proposal)
            use pm_matrixChol, only: setMatChol, lowDia, transHerm
            use pm_matrixTrace, only: getMatMulTraceLog
            use pm_ellipsoid, only: getLogVolUnitBall ! for uniform proposal.
            use pm_arrayRebind, only: setRebound
            type(err_type), intent(inout) :: err
            type(spec_type), intent(in) :: spec
            type(proposal_type) :: proposal
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@constructProposal()"
            integer(IK) :: info
            !call setNAN(proposal%logSqrtDetOld)
            ndim = spec%ndim%val
            ndimp1 = ndim + 1_IK
            proposal%format = SKC_"("//spec%ndim%str//SKC_"(g0,:,', '),"//new_line(SKC_"a")//SKC_")"
            proposal%meanOld = spec%proposalStart%val
            proposal%sampleSizeOld = 1_IK ! This initialization is crucial important for proposal adaptation.
            proposal%scaleFactorSq%running = 1._RKC
            proposal%scaleFactorSq%default = spec%proposalScaleFactor%val**2
            proposal%delRejEnabled = 0_IK < spec%proposalDelayedRejectionCount%val
            !proposal%delayedRejection%scaleFactor%val = real(spec%proposalDelayedRejectionScaleFactor%val, RKC)
            !proposal%delayedRejection%scaleFactor%log = log(proposal%delayedRejection%scaleFactor%val)
            !proposal%domainCubeLimitLower = real(spec%domainCubeLimitLower%val, RKC)
            !proposal%domainCubeLimitUpper = real(spec%domainCubeLimitUpper%val, RKC)
            !if (spec%domain%isFinite) proposal%domainBallCovMatInv = real(spec%domainBallCovMat%inv, RKC)
            !proposal%domainBallCenter = real(spec%domainBallCenter%val, RKC)
            ! setup covariance matrix
            SET_ParaDISE(call setRebound(proposal%logSqrtDetInvCov, 0_IK, spec%proposalDelayedRejectionCount%val))
            SET_ParaDISE(call setRebound(proposal%invCov, [1_IK, 1_IK, 0_IK], [ndim, ndimp1, spec%proposalDelayedRejectionCount%val]))
            call setRebound(proposal%covLowChoUpp, [1_IK, 1_IK, 0_IK], [ndim, ndimp1, spec%proposalDelayedRejectionCount%val])
            SET_CAF(if (allocated(comv_covLowChoUpp)) deallocate(comv_covLowChoUpp))
            SET_CAF(allocate(comv_covLowChoUpp(ndim, 1 : ndimp1)[*]))
            proposal%covLowChoUpp(1 : ndim, 1 : ndim, 0) = spec%proposalCovMat%val * proposal%scaleFactorSq%default ! Scale the initial covariance matrix
            call setMatChol(proposal%covLowChoUpp(:, 1 : ndim, 0), lowDia, info, proposal%covLowChoUpp(:, 2 : ndimp1, 0), transHerm) ! The `:` instead of `1 : ndim` avoids temporary array creation.
            err%occurred = info /= 0_IK
            if (err%occurred) then
                err%msg = "Internal erorr: Cholesky factorization of proposalCovMat failed at diagonal #"//getStr(info)//SKC_": "//getStr(transpose(proposal%covLowChoUpp(:, 1 : ndim, 0)), proposal%format)
                return
            end if
            proposal%logSqrtDetOld = .5_RKC * getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : ndimp1, 0))
            ! Scale the higher-stage delayed-rejection Cholesky Upper matrices.
            if (proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal)
            SET_ParaDISE(proposal%negLogVolUnitBall = -getLogVolUnitBall(real(ndim, RKC))) ! only for uniform proposal?
            SET_ParaDISE(call setProposalInvCov(spec, proposal))
            ! read/write the first entry of the restart file
            if (spec%image%is%leader) then
                block
                    real(RKC) :: meanAccRateSinceStart
                    if (spec%run%is%new) then
                        !if (.not. spec%outputRestartFileFormat%isBinary) then
                        !    write(spec%restartFile%unit) "ndim"
                        !    write(spec%restartFile%unit)  ndim
                        !end if
                        call writeRestart(spec, proposal, meanAccRateSinceStart = 1._RKC)
                    else
                        !if (.not. spec%outputRestartFileFormat%isBinary) then
                        !    read(spec%restartFile%unit)
                        !    read(spec%restartFile%unit)
                        !end if
                        call readRestart(spec, meanAccRateSinceStart)
                    end if
                end block
            end if
        end function constructProposal

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !>  Update the Cholesky factorization of the higher-stage delayed rejection covariance matrices.
        !>
        !>  \todo
        !>  The performance of this update could be improved by only updating the higher-stage covariance, only when needed.<br>
        !>  However, the gain will be likely minimal, especially in low-dimensions.<br>
        pure subroutine setProposalDelRejChol(spec, proposal)
            type(proposal_type), intent(inout) :: proposal
            type(spec_type), intent(in) :: spec
            integer(IK) :: idim, ndim, istage
            ndim = size(proposal%covLowChoUpp, 1, IK)
            do istage = 1, spec%proposalDelayedRejectionCount%val
                ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                !do concurrent(idim = 1 : ndim)
                do idim = 1, ndim
                    proposal%covLowChoUpp(1 : idim, idim + 1, istage) = proposal%covLowChoUpp(1 : idim, idim + 1, istage - 1) * spec%proposalDelayedRejectionScaleFactor%val(istage)
                end do
            end do
            ! There is no need to check for positive-definiteness of the covLowChoUpp, it is already checked on the first image.
        end subroutine setProposalDelRejChol

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine readRestart(spec, meanAccRateSinceStart)
            real(RKC), intent(out) :: meanAccRateSinceStart
            type(spec_type), intent(in) :: spec
            integer(IK) :: idim
            if (spec%outputRestartFileFormat%isBinary) then
                read(spec%restartFile%unit) meanAccRateSinceStart
            else
                read(spec%restartFile%unit, *)
                read(spec%restartFile%unit, *) meanAccRateSinceStart
                do idim = 1, 8 + ndim * (ndim + 3) / 2; read(spec%restartFile%unit, *); end do
            end if
        end subroutine readRestart

        subroutine writeRestart(spec, proposal, meanAccRateSinceStart)
            real(RKC), intent(in) :: meanAccRateSinceStart
            type(proposal_type), intent(in) :: proposal
            type(spec_type), intent(in) :: spec
            integer(IK) :: idim, jdim
            if (spec%outputRestartFileFormat%isBinary) then
                write(spec%restartFile%unit) meanAccRateSinceStart
            else
                ! The extra components are all informational for the end user, otherwise not needed for a restart.
                write(spec%restartFile%unit, spec%restartFile%format) "meanAcceptanceRateSinceStart" & ! LCOV_EXCL_LINE
                                                                    , meanAccRateSinceStart & ! LCOV_EXCL_LINE
                                                                    , "numFuncCall" & ! LCOV_EXCL_LINE
                                                                    , proposal%sampleSizeOld & ! LCOV_EXCL_LINE
                                                                    , "proposalAdaptiveScaleFactorSquared" & ! LCOV_EXCL_LINE
                                                                    , proposal%scaleFactorSq%running * proposal%scaleFactorSq%default & ! LCOV_EXCL_LINE
                                                                    , "proposalLogVolume" & ! LCOV_EXCL_LINE
                                                                    , proposal%logSqrtDetOld & ! LCOV_EXCL_LINE
                                                                    , "proposalMean" & ! LCOV_EXCL_LINE
                                                                    , proposal%meanOld & ! LCOV_EXCL_LINE
                                                                    , "proposalCov" & ! LCOV_EXCL_LINE
                                                                    , ((proposal%covLowChoUpp(idim, jdim, 0), idim = jdim, ndim), jdim = 1, ndim)
                                                                    ! Only the lower-triangle of the covariance matrix with no delayed rejection.
            end if
            flush(spec%restartFile%unit)
        end subroutine writeRestart

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Return a newly sampled state to the kernel routine for acceptance check.
        ! ATTN: The colon index in place of 1 : ndim below avoids the temporary array creation.
        ! The purity breaks only because of the warn() call.
        subroutine setProposalStateNew(spec, proposal, delayedRejectionStage, stateOld, stateNew, rng, err)
            use pm_ellipsoid, only: isMemberEll
            use pm_distUnifEll, only: setUnifEllRand
            use pm_distMultiNorm, only: setMultiNormRand, uppDia
            use pm_distUnif, only: xoshiro256ssw_type
            type(proposal_type), intent(inout) :: proposal
            type(xoshiro256ssw_type), intent(inout) :: rng
            real(RKC), intent(in), contiguous :: stateOld(:)
            real(RKC), intent(out), contiguous :: stateNew(:)
            integer(IK), intent(in) :: delayedRejectionStage
            type(spec_type), intent(inout) :: spec
            type(err_type), intent(inout) :: err
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@setProposalStateNew()"
            character(*, SK), parameter :: WARNING_MSG = PROCEDURE_NAME//SKC_"Number of proposed states out of the objective function domain without any acceptance: "
            character(*, SK), parameter :: FATAL_PREFIX = PROCEDURE_NAME//SKC_"Number of proposed states out of the objective function domain without any acceptance has reach the maximum specified value (domainErrCountMax = "
            character(*, SK), parameter :: FATAL_SUFFIX = SKC_"). This can occur if the domain of the density function is incorrectly specified or if the density function definition implementation or definition is flawed."
            integer(IK) :: domainCheckCounter
            domainCheckCounter = 1_IK
            loopBoundaryCheck: do ! Check for the support Region consistency:
                ! \todo
                ! There is an extremely low but non-zero likelihood that the resulting `stateNew` would be the same as `stateOld`.
                ! This could later become problematic because the same state is accepted as new state, potentially corrupting
                ! the proposal covariance computation based on the newly-collected samples.
                if (spec%proposal%is%normal) then
                    call setMultiNormRand(rng, stateNew, stateOld, proposal%covLowChoUpp(:, 2 : ndimp1, delayedRejectionStage), uppDia)
                elseif (spec%proposal%is%uniform) then
                    call setUnifEllRand(rng, stateNew, stateOld, proposal%covLowChoUpp(:, 2 : ndimp1, delayedRejectionStage), uppDia)
                else
                    error stop "Internal error occurred: proposal distribution unrecognized."
                end if
                ! check if proposal is in the domain.
                if (spec%domain%isFinite) then
                    if (spec%domain%isCube) then
                        if (all(spec%domainCubeLimitLower%val <= stateNew .and. stateNew <= spec%domainCubeLimitUpper%val)) return
                    elseif (spec%domain%isBall) then
                        if (isMemberEll(spec%domainBallCovMat%inv, spec%domainBallCenter%val, stateNew)) return
                    end if
                else
                    return
                end if
                if (domainCheckCounter == spec%domainErrCount%val) call spec%disp%warn%show(WARNING_MSG//getStr(domainCheckCounter))
                if (domainCheckCounter == spec%domainErrCountMax%val) exit loopBoundaryCheck
                domainCheckCounter = domainCheckCounter + 1_IK
            end do loopBoundaryCheck
            err%msg = FATAL_PREFIX//spec%domainErrCountMax%str//FATAL_SUFFIX
            err%occurred = .true._LK
        end subroutine setProposalStateNew

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !>  Adapt the proposal distribution based on the specified weighted sampleState.
        !>
        !>  \param[in]      sampleState             :   The two dimensional `(ndim, nsam)` array of accepted states.
        !>  \param[in]      sampleWeight            :   The one dimensional integer array of the weights of the sampled states in the sampleState.
        !>  \param[in]      adaptationIsGreedy      :   The logical flag indicating whether the adaptation is in greedy mode..
        !>  \param[inout]   meanAccRateSinceStart   :   The mean acceptance rate since the start of the sampling.
        !>  \param[out]     adaptationSucceeded     :   The output logical flag indicating whether the sampling succeeded.
        !>  \param[out]     adaptationMeasure       :   The output real number in the range `[0,1]` indicating the amount of adaptation,
        !>                                              with zero indicating no adaptation and one indicating extreme adaptation to the extent
        !>                                              that the new adapted proposal distribution is completely different from the previous proposal.
        !>  \warning
        !>  This routine must be exclusively called by the leader images.
        !>
        !>  \note
        !>  Other than the call to `warn%show()`, this procedure is `pure`.<br>
        subroutine setProposalAdapted(spec, proposal, sampleState, sampleWeight, adaptationIsGreedy, meanAccRateSinceStart, adaptationSucceeded, adaptationMeasure, err)
            use pm_sampleMean, only: setMean
            use pm_matrixTrace, only: getMatMulTraceLog
            use pm_sampleCov, only: setCov, setCovMeanMerged
            use pm_matrixChol, only: setMatChol, lowDia, transHerm
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@setProposalAdapted()"
            type(spec_type), intent(inout) :: spec
            type(proposal_type), intent(inout) :: proposal
            real(RKC), intent(in), contiguous :: sampleState(:,:)
            integer(IK), intent(in), contiguous :: sampleWeight(:)
            logical(LK), intent(in) :: adaptationIsGreedy
            real(RKC), intent(inout) :: meanAccRateSinceStart ! is intent(out) in restart mode, intent(in) in fresh mode.
            logical(LK), intent(out) :: adaptationSucceeded
            real(RKC), intent(out) :: adaptationMeasure
            type(err_type), intent(inout) :: err
            real(RKC) :: sampleStateMeanNew(size(sampleState, 1, IK))
            real(RKC) :: sampleStateMeanCurrent(size(sampleState, 1, IK))
            real(RKC) :: sampleStateCovCholOld(size(sampleState, 1, IK), size(sampleState, 1, IK) + 1)
            real(RKC) :: sampleStateCovLowCurrent(size(sampleState, 1, IK), size(sampleState, 1, IK) + 1)
            real(RKC) :: logSqrtDetNew, logSqrtDetSum, adaptiveScaleFactor, fracA
            logical(LK) :: adaptationMeasureComputationNeeded ! only used to avoid redundant affinity computation, if no update occurs.
            logical(LK) :: scalingNeeded!, singularityOccurred
            integer(IK) :: sampleSizeOld, sampleSizeCurrent
            integer(IK):: idim, jdim, ndim, nsam, info
            integer(IK), parameter :: dim = 2_IK
            scalingNeeded = .false._LK
            ndim = size(sampleState, 1, IK)
            nsam = size(sampleState, 2, IK)
            sampleSizeOld = proposal%sampleSizeOld ! this is kept only for restoration of proposal%sampleSizeOld, if needed.
            if (spec%image%is%leader .and. spec%run%is%dry) call readRestart(spec, meanAccRateSinceStart)
            ! First if there are less than ndimp1 points for new covariance computation, then just scale the covariance and return.
            blockSufficientSampleSizeCheck: if (ndim < nsam) then
                ! get the upper covariance matrix and mean of the new sample
                if (adaptationIsGreedy) then
                    sampleSizeCurrent = nsam
                    call setMean(sampleStateMeanCurrent, sampleState, dim)
                    call setCov(sampleStateCovLowCurrent(:, 1 : ndim), lowDia, sampleStateMeanCurrent, sampleState, dim)
                else
                    call setMean(sampleStateMeanCurrent, sampleState, dim, sampleWeight, sampleSizeCurrent)
                    call setCov(sampleStateCovLowCurrent(:, 1 : ndim), lowDia, sampleStateMeanCurrent, sampleState, dim, sampleWeight, sampleSizeCurrent)
                end if
                ! Combine old and new covariance matrices if both exist.
                sampleStateCovCholOld = proposal%covLowChoUpp(:, :, 0) ! This will be used to recover the old covariance in case of update failure, and to compute the adaptation measure.
                blockMergeCovMat: if (proposal%sampleSizeOld == 1_IK) then
                    ! There is no prior old Covariance matrix to combine with the new one from the new sampleState
                    proposal%meanOld(1 : ndim) = sampleStateMeanCurrent
                    proposal%sampleSizeOld = sampleSizeCurrent
                    ! Copy and then scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor.
                    ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                    !do concurrent(jdim = 1 : ndim)
                    do jdim = 1, ndim
                        proposal%covLowChoUpp(jdim : ndim, jdim, 0) = sampleStateCovLowCurrent(jdim : ndim, jdim) * proposal%scaleFactorSq%default
                    end do
                else blockMergeCovMat
                    ! First scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor.
                    ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                    !do concurrent(jdim = 1 : ndim)
                    do jdim = 1, ndim
                        sampleStateCovLowCurrent(jdim : ndim, jdim) = sampleStateCovLowCurrent(jdim : ndim, jdim) * proposal%scaleFactorSq%default
                    end do
                    ! Now combine it with the old covariance matrix.
                    ! Do not set the full boundaries' range `(1 : ndim)` for `proposal%covLowChoUpp` in the following subroutine call.
                    ! Setting the boundaries forces the compiler to generate a temporary array.
                    fracA = real(proposal%sampleSizeOld, RKC) / real(proposal%sampleSizeOld + sampleSizeCurrent, RKC)
                    call setCovMeanMerged(proposal%covLowChoUpp(:, 1 : ndim, 0), sampleStateMeanNew, sampleStateCovLowCurrent(:, 1 : ndim), sampleStateMeanCurrent, sampleStateCovCholOld(:, 1 : ndim), proposal%meanOld, fracA, lowDia)
                    proposal%meanOld(1 : ndim) = sampleStateMeanNew
                end if blockMergeCovMat
                ! Now get the Cholesky factorization.
                ! WARNING: Do not set the full boundaries' range `(1 : ndim)` for the first index of `proposal%covLowChoUpp` in the following subroutine call.
                ! WARNING: Setting the boundaries forces the compiler to generate a temporary array.
                call setMatChol(proposal%covLowChoUpp(:, 1 : ndim, 0), lowDia, info, proposal%covLowChoUpp(:, 2 : ndimp1, 0), transHerm) ! The `:` instead of `1 : ndim` avoids temporary array creation.
                blockPosDefCheck: if (info == 0_IK) then
                    !singularityOccurred = .false._LK
                    adaptationSucceeded = .true._LK
                    adaptationMeasureComputationNeeded = .true._LK
                    proposal%sampleSizeOld = proposal%sampleSizeOld + sampleSizeCurrent
                else blockPosDefCheck
                    adaptationMeasure = 0._RKC
                    !singularityOccurred = .true._LK
                    adaptationSucceeded = .false._LK
                    adaptationMeasureComputationNeeded = .false._LK
                    proposal%sampleSizeOld = sampleSizeOld
                    ! It may be a good idea to add a warning message printed out here for the singularity occurrence.
                    call spec%disp%warn%show(SKC_"Singularity occurred while adapting the proposal distribution covariance matrix.", tmsize = 0_IK, bmsize = 0_IK)
                    ! Recover the old covariance matrix + Cholesky factor.
                    proposal%covLowChoUpp(:, :, 0) = sampleStateCovCholOld
                end if blockPosDefCheck
                !! perform global adaptive scaling is requested
                !if (spec%targetAcceptanceRate%enabled) then
                !    if (meanAccRateSinceStart < spec%targetAcceptanceRate%aim) then
                !        proposal%scaleFactorSq%running = mc_maxScaleFactorSq**((spec%targetAcceptanceRate%aim-meanAccRateSinceStart)/spec%targetAcceptanceRate%aim)
                !    else
                !        proposal%scaleFactorSq%running = mc_maxScaleFactorSq**((spec%targetAcceptanceRate%aim-meanAccRateSinceStart)/(1._RKC-spec%targetAcceptanceRate%aim))
                !    end if
                !    proposal%scaleFactorSq%running = proposal%scaleFactorSq%running * (meanAccRateSinceStart/spec%targetAcceptanceRate%aim)**spec%ndim%inv
                !end if
                if (spec%targetAcceptanceRate%enabled) scalingNeeded = .true._LK
            else blockSufficientSampleSizeCheck
                ! Singularity has occurred. If the first covariance merging has not occurred yet, set the scaling factor appropriately to shrink the covariance matrix.
                adaptationSucceeded = .false._LK
                if (proposal%sampleSizeOld == 1_IK .or. spec%targetAcceptanceRate%enabled) then
                    scalingNeeded = .true._LK
                    adaptationMeasureComputationNeeded = .true._LK
                    ! Save the old covariance matrix for the computation of the adaptation measure.
                    ! \todo Could these copies be somehow removed and optimized away? perhaps unimportant if sample size is sufficient most of the times.
                    do jdim = 1, ndim
                        do idim = jdim, ndim
                            sampleStateCovCholOld(idim, jdim) = proposal%covLowChoUpp(idim, jdim, 0)
                        end do
                    end do
                else
                    adaptationMeasure = 0._RKC
                    adaptationMeasureComputationNeeded = .false._LK
                end if
            end if blockSufficientSampleSizeCheck
            ! Adjust the scale of the covariance matrix and the Cholesky factor, if needed.
            !scalingNeeded = .true._LK
            if (scalingNeeded) then
                if (meanAccRateSinceStart < spec%targetAcceptanceRate%val(1) .or. spec%targetAcceptanceRate%val(2) < meanAccRateSinceStart) then
                    proposal%scaleFactorSq%running = (meanAccRateSinceStart / spec%targetAcceptanceRate%aim)**spec%ndim%inv
                    !block
                    !    use pm_distUnif, only: getUnifRand
                    !    integer, save :: counter = 0_IK
                    !    counter = counter - 1
                    !    proposal%scaleFactorSq%running = proposal%scaleFactorSq%running * exp(-counter*getUnifRand(-1._RKC, 1._RKC)/1.e4_RK)
                    !    !use pm_statistics, only: getUnifRand
                    !    !proposal%scaleFactorSq%running = proposal%scaleFactorSq%running * exp(real(getUnifRand(-1_IK, 1_IK), RK))
                    !    !write(*,*) counter, proposal%scaleFactorSq%running
                    !end block
                    adaptiveScaleFactor = sqrt(proposal%scaleFactorSq%running)
                    proposal%covLowChoUpp(1 : ndim, 1, 0) = proposal%covLowChoUpp(1 : ndim, 1, 0) * proposal%scaleFactorSq%running ! Update first column of covariance matrix.
                    ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                    !do concurrent(jdim = 2 : ndim)
                    do jdim = 2, ndim
                        proposal%covLowChoUpp(1 : jdim - 1, jdim, 0) = proposal%covLowChoUpp(1 : jdim, jdim + 1, 0) * adaptiveScaleFactor ! Update the Cholesky factor.
                        proposal%covLowChoUpp(jdim : ndim, jdim, 0) = proposal%covLowChoUpp(jdim : ndim, jdim, 0) * proposal%scaleFactorSq%running ! Update covariance matrix.
                    end do
                    proposal%covLowChoUpp(1 : ndim, ndimp1, 0) = proposal%covLowChoUpp(1 : ndim, ndimp1, 0) * adaptiveScaleFactor ! Update last column of Cholesky factor.
                end if
            end if
            ! Compute the adaptivity only if any updates has occurred.
            blockAdaptationMeasureComputation: if (adaptationMeasureComputationNeeded) then
                logSqrtDetNew = getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : ndimp1, 0))
                ! Use the universal upper bound.
                !block
                !use pm_io
                !call disp%show(sampleStateCovLowCurrent)
                !end block
                ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                !do concurrent(jdim = 1 : ndim)
                do jdim = 1, ndim
                    sampleStateCovLowCurrent(jdim : ndim, jdim) = 0.5_RKC * (proposal%covLowChoUpp(jdim : ndim, jdim, 0) + sampleStateCovCholOld(jdim : ndim, jdim)) ! dummy
                end do
                call setMatChol(sampleStateCovLowCurrent(:, 1 : ndim), lowDia, info, sampleStateCovLowCurrent(:, 2 : ndimp1), transHerm)
                if (info == 0_IK) then !singularityOccurred = info /= 0_IK
                    logSqrtDetSum = getMatMulTraceLog(sampleStateCovLowCurrent(:, 2 : ndimp1))
                else ! singularityOccurred ! LCOV_EXCL_START
                    err%occurred = .true._LK
                    err%msg = NL1//PROCEDURE_NAME//SKC_"@line::"//getStr(__LINE__)//SKC_": "//&
                    SKC_"Singular lower covariance matrix detected at column ("//getStr(info)//SKC_") while computing the adaptationMeasure measure:"//NL2
                    do info = 1, size(sampleStateCovLowCurrent, 1, IK)
                        err%msg = err%msg//getStr(sampleStateCovLowCurrent(info,:))//NL1
                    end do
                    err%msg = err%msg//NL1//&
                    SKC_": Error occurred while computing the Cholesky factorization of a matrix needed for the computation of the adaptation measure. &
                    &Such error is highly unusual, and requires an in depth investigation of the case. &
                    &It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations. &
                    &In such case, consider increasing the value of the input variable proposalAdaptationPeriod. &
                    &It may also be that your input objective function has been incorrectly implemented. &
                    &Also, ensure a correct value of `ndim` to the ParaMonte sampler routine, representing the domain size. &
                    &Otherwise, restarting the simulation might resolve the error. If the error persists, please inform the developers at: "//NL2//PARAMONTE_WEB_ISSUES
                    return
                end if ! LCOV_EXCL_STOP
                adaptationMeasure = 1._RKC - exp(0.5_RKC * (proposal%logSqrtDetOld + logSqrtDetNew) - logSqrtDetSum) ! totalVariationUpperBound
                if (0._RKC < adaptationMeasure) then
                    adaptationMeasure = sqrt(adaptationMeasure) ! totalVariationUpperBound
                elseif (adaptationMeasure < 0._RKC) then
                    spec%msg = PROCEDURE_NAME//SKC_": Non-positive adaptation measure detected" ! LCOV_EXCL_LINE
                    if (abs(adaptationMeasure) < sqrt(epsilon(0._RKC))) spec%msg = spec%msg//SKC_" (possibly due to round-off error)" ! LCOV_EXCL_LINE
                    call spec%disp%warn%show(spec%msg//SKC_": "//getStr(adaptationMeasure)) ! LCOV_EXCL_LINE
                    adaptationMeasure = 0._RKC ! LCOV_EXCL_LINE
                end if
                proposal%logSqrtDetOld = logSqrtDetNew
                !block
                !integer, save :: counter = 0
                !counter = counter + 1
                !!if (counter==1) then
                !if (adaptationMeasure>1._RKC) then
                !write(*,*)
                !write(*,*) proposal%logSqrtDetOld
                !write(*,*) logSqrtDetNew
                !write(*,*) logSqrtDetSum
                !write(*,*) proposal%logSqrtDetOld + logSqrtDetNew - 2_IK * logSqrtDetSum
                !write(*,*) exp( proposal%logSqrtDetOld + logSqrtDetNew - 2_IK * logSqrtDetSum )
                !write(*,*)
                !end if
                !end block
                ! update the higher-stage delayed-rejection Cholesky Upper matrices
                if (proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal)
                SET_ParaDISE(call setProposalInvCov(spec, proposal))
            end if blockAdaptationMeasureComputation
            if (spec%image%is%leader .and. spec%run%is%new) call writeRestart(spec, proposal, meanAccRateSinceStart)
        end subroutine setProposalAdapted

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  RETURN_IF_FAILED
#undef  SET_ParaDISE

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef SET_CAF_MPI
#undef SET_CAF
#undef SHARED
