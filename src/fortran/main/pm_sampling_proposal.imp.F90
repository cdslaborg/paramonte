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
err%msg = PROCEDURE_NAME//getLine(LINE)//SKG_": "//MSG; \
call spec%disp%stop%show(err%msg); \
return; \
end if;
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     ParaDISE_ENABLED || ParaDRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_err, only: getLine
        use pm_err, only: err_type
        use pm_val2str, only: getStr
        use pm_kind, only: SKG => SK, SK, IK, LK
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
        character(*,SKG), parameter :: MODULE_NAME = SK_"@pm_sampling_proposal_dise"
#elif   ParaDRAM_ENABLED
        use pm_sampling_dram, only: NL1, NL2
        use pm_sampling_dram, only: spec_type => specdram_type
        use pm_sampling_dram, only: stat_type => statdram_type
        use pm_sampling_scio, only: getErrChainRead => getErrChainReadDRAM
        use pm_sampling_scio, only: getErrChainWrite => getErrChainWriteDRAM
        use pm_sampling_scio, only: chainFileColName => chainFileColNameDRAM
        use pm_sampling_scio, only: isFailedChainResize => isFailedChainResizeDRAM
        implicit none
        character(*,SKG), parameter :: MODULE_NAME = SK_"@pm_sampling_proposal_dram"
#endif
#if     CAF_ENABLED
        real(RKG), save     , allocatable   :: comv_covLowChoUpp(:,:)[:] ! lower covariance, upper Cholesky.
#endif
        type :: scaleSq_type
            real(RKG) :: default
            real(RKG) :: running
        end type
        type :: proposal_type
            ! the following are made components for the sake of thread-safe save attribute and the restart file generation.
            character(:,SKG), allocatable   :: format
            type(scaleSq_type)              :: scaleSq
            logical(LK)                     :: delRejEnabled
            integer(IK)                     :: sampleSizeOld
            real(RKG)                       :: logSqrtDetOld
            real(RKG)       , allocatable   :: covLowChoUpp(:,:,:) ! The covariance Matrix of the proposal distribution. Last index belongs to delayed rejection.
            real(RKG)       , allocatable   :: meanOld(:)
#if         ParaDISE_ENABLED
            real(RKG)                       :: negLogVolUnitBall
            real(RKG)       , allocatable   :: logSqrtDetInvCov(:) ! (0 : spec%proposalDelayedRejectionCount%val)
            real(RKG)       , allocatable   :: invCov(:,:,:)
#define     SET_ParaDISE(X)X
#else
#define     SET_ParaDISE(X)
#endif
        end type
        interface proposal_type
            procedure :: proposal_typer
        end interface

contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine killMeAlreadyCMake1_RK5(); use pm_sampling_scio_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake1_RK4(); use pm_sampling_scio_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake1_RK3(); use pm_sampling_scio_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake1_RK2(); use pm_sampling_scio_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake1_RK1(); use pm_sampling_scio_RK1, only: RKG; end subroutine
#if     ParaDISE_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_dise_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_dise_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_dise_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_dise_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_dise_RK1, only: RKG; end subroutine
#elif   ParaDRAM_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_dram_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_dram_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_dram_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_dram_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_dram_RK1, only: RKG; end subroutine
#elif   ParaNest_ENABLED
        subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_nest_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_nest_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_nest_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_nest_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_nest_RK1, only: RKG; end subroutine
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
                comv_covLowChoUpp(1 : spec%ndim%val, 1 : spec%ndim%vp1) = proposal%covLowChoUpp(1: spec%ndim%val, 1 : spec%ndim%vp1, 0)
            else
                proposal%covLowChoUpp(1: spec%ndim%val, 1 : spec%ndim%vp1, 0) = comv_covLowChoUpp(1 : spec%ndim%val, 1 : spec%ndim%vp1)[1]
                if (proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal) ! update the higher-stage delayed-rejection Cholesky Upper matrices.
                SET_ParaDISE(call setProposalInvCov(proposal))
            end if
#elif       MPI_ENABLED
            use mpi !mpi_f08, only: mpi_bcast, mpi_double_precision, mpi_comm_world ! LCOV_EXCL_LINE
            use, intrinsic :: iso_c_binding, only: c_sizeof
            integer, parameter :: NBRK = c_sizeof(0._RKG) ! Number of bytes in the real kind.
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
            call setMatInv(proposal%invCov(:, 1 : spec%ndim%val, istage), proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, istage), choUpp)
            proposal%logSqrtDetInvCov(istage) = -proposal%logSqrtDetOld
            do istage = 0, spec%proposalDelayedRejectionCount%val
                call setMatInv(proposal%invCov(:, 1 : spec%ndim%val, istage), proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, istage), choUpp)
                !proposal%logSqrtDetInvCov(istage) = .5_RKG * getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, istage))
                proposal%logSqrtDetInvCov(istage) = proposal%logSqrtDetInvCov(istage) - spec%ndim%val * log(spec%proposalDelayedRejectionScale%val(istage))
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
            real(RKG), intent(in), contiguous :: origin(:), destin(:)
            character(*,SKG), parameter :: MSG = SKG_"@getProposalJumpLogProb()/: Internal library error; The specified uniform proposal jump is not within the proposal domain. "
            character(*,SKG), parameter :: ERRMSG = MSG//SKG_"ndim, delayedRejectionStage, proposal%invCov(:,1 : ndim, delayedRejectionStage), origin, destin = "
            real(RKG) :: logPDF
            if (spec%proposal%is%normal) then
                logPDF = getMultiNormLogPDF(destin, origin, proposal%invCov(:,1 : spec%ndim%val, delayedRejectionStage), logPDFNF = getMultiNormLogPDFNF(spec%ndim%val, proposal%logSqrtDetInvCov(delayedRejectionStage)))
            elseif (spec%proposal%is%uniform) then
                logPDF = proposal%negLogVolUnitBall + proposal%logSqrtDetInvCov(delayedRejectionStage)
                CHECK_ASSERTION(__LINE__, isMemberEll(proposal%invCov(:,1 : spec%ndim%val, delayedRejectionStage), origin, destin), \
                ERRMSG//getStr([spec%ndim%val, delayedRejectionStage])//SKG_", "//getStr([proposal%invCov(:,1 : spec%ndim%val, delayedRejectionStage), origin, destin]))
            end if
        end function getProposalJumpLogProb
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function proposal_typer(spec, err) result(proposal)
            use pm_matrixChol, only: setMatChol, lowDia, transHerm
            use pm_matrixTrace, only: getMatMulTraceLog
            use pm_ellipsoid, only: getLogVolUnitBall ! for uniform proposal.
            use pm_arrayRebind, only: setRebound
            type(err_type), intent(inout) :: err
            type(spec_type), intent(in) :: spec
            type(proposal_type) :: proposal
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@proposal_typer()"
            integer(IK) :: info
            !call setNAN(proposal%logSqrtDetOld)
            proposal%format = SKG_"("//spec%ndim%str//SKG_"(g0,:,', '),"//new_line(SKG_"a")//SKG_")"
            proposal%meanOld = spec%proposalStart%val
            proposal%sampleSizeOld = 1_IK ! This initialization is crucial important for proposal adaptation.
            proposal%scaleSq%running = 1._RKG
            proposal%scaleSq%default = spec%proposalScale%val**2
            proposal%delRejEnabled = 0_IK < spec%proposalDelayedRejectionCount%val
            !proposal%delayedRejection%scale%val = real(spec%proposalDelayedRejectionScale%val, RKG)
            !proposal%delayedRejection%scale%log = log(proposal%delayedRejection%scale%val)
            !proposal%domainCubeLimitLower = real(spec%domainCubeLimitLower%val, RKG)
            !proposal%domainCubeLimitUpper = real(spec%domainCubeLimitUpper%val, RKG)
            !if (spec%domain%isFinite) proposal%domainBallCovInv = real(spec%domainBallCov%inv, RKG)
            !proposal%domainBallAvg = real(spec%domainBallAvg%val, RKG)
            ! setup covariance matrix
            SET_ParaDISE(call setRebound(proposal%logSqrtDetInvCov, 0_IK, spec%proposalDelayedRejectionCount%val))
            SET_ParaDISE(call setRebound(proposal%invCov, [1_IK, 1_IK, 0_IK], [spec%ndim%val, spec%ndim%vp1, spec%proposalDelayedRejectionCount%val]))
            call setRebound(proposal%covLowChoUpp, [1_IK, 1_IK, 0_IK], [spec%ndim%val, spec%ndim%vp1, spec%proposalDelayedRejectionCount%val])
            SET_CAF(if (allocated(comv_covLowChoUpp)) deallocate(comv_covLowChoUpp))
            SET_CAF(allocate(comv_covLowChoUpp(spec%ndim%val, 1 : spec%ndim%vp1)[*]))
            proposal%covLowChoUpp(1 : spec%ndim%val, 1 : spec%ndim%val, 0) = spec%proposalCov%val * proposal%scaleSq%default ! Scale the initial covariance matrix
            ! The `:` instead of `1 : spec%ndim%val` avoids temporary array creation.
            call setMatChol(proposal%covLowChoUpp(:, 1 : spec%ndim%val, 0), lowDia, info, proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, 0), transHerm)
            err%occurred = info /= 0_IK
            if (err%occurred) then
                err%msg = "Internal erorr: Cholesky factorization of proposalCov failed at diagonal #"// & ! LCOV_EXCL_LINE
                getStr(info)//SKG_": "//getStr(transpose(proposal%covLowChoUpp(:, 1 : spec%ndim%val, 0)), proposal%format)
                return
            end if
            proposal%logSqrtDetOld = .5_RKG * getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, 0))
            ! Scale the higher-stage delayed-rejection Cholesky Upper matrices.
            if (proposal%delRejEnabled) call setProposalDelRejChol(spec, proposal)
            SET_ParaDISE(proposal%negLogVolUnitBall = -getLogVolUnitBall(real(spec%ndim%val, RKG))) ! only for uniform proposal?
            SET_ParaDISE(call setProposalInvCov(spec, proposal))
            ! read/write the first entry of the restart file
            if (spec%image%is%leader) then
                block
                    integer(IK) :: idim
                    real(RKG) :: meanAccRateSinceStart
                    if (spec%run%is%new) then
                        if (.not. spec%outputRestartFileFormat%isBinary) then
                            write(spec%restartFile%unit, spec%restartFile%format) "ndim"
                            write(spec%restartFile%unit, spec%restartFile%format)  spec%ndim%val
                            write(spec%restartFile%unit, spec%restartFile%format)  "domainAxisName"
                            do idim = 1, spec%ndim%val
                                write(spec%restartFile%unit, spec%restartFile%format) trim(adjustl(spec%domainAxisName%val(idim)))
                            end do
                        end if
                        call writeRestart(spec, meanAccRateSinceStart = 1._RKG, proposal = proposal)
                    else
                        if (.not. spec%outputRestartFileFormat%isBinary) then
                            do idim = 1, spec%ndim%val + 3
                                read(spec%restartFile%unit)
                            end do
                        end if
                        call readRestart(spec, meanAccRateSinceStart)
                    end if
                end block
            end if
        end function proposal_typer

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
            integer(IK) :: idim, istage
            do istage = 1, spec%proposalDelayedRejectionCount%val
                ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                !do concurrent(idim = 1 : spec%ndim%val)
                do idim = 1, spec%ndim%val
                    proposal%covLowChoUpp(1 : idim, idim + 1, istage) = proposal%covLowChoUpp(1 : idim, idim + 1, istage - 1) * spec%proposalDelayedRejectionScale%val(istage)
                end do
            end do
            ! There is no need to check for positive-definiteness of the covLowChoUpp, it is already checked on the first image.
        end subroutine setProposalDelRejChol

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine readRestart(spec, meanAccRateSinceStart, nskip)
            real(RKG), intent(out) :: meanAccRateSinceStart
            integer(IK), intent(in), optional :: nskip
            type(spec_type), intent(in) :: spec
            integer(IK) :: nskip_def, idim
            if (spec%outputRestartFileFormat%isBinary) then
                read(spec%restartFile%unit) meanAccRateSinceStart
            else
                read(spec%restartFile%unit, *)
                read(spec%restartFile%unit, *) meanAccRateSinceStart
                nskip_def = 8 + spec%ndim%val * (spec%ndim%val + 3) / 2
                if (present(nskip)) nskip_def = nskip
                !do idim = 1, 8 + spec%ndim%val * (spec%ndim%val + 3) / 2 !nskip_def
                do idim = 1, nskip_def
                    read(spec%restartFile%unit, *)
                end do
            end if
        end subroutine readRestart

        subroutine writeRestart(spec, meanAccRateSinceStart, proposal)
            type(proposal_type), intent(in), optional :: proposal
            real(RKG), intent(in) :: meanAccRateSinceStart
            type(spec_type), intent(in) :: spec
            integer(IK) :: idim, jdim
            if (spec%outputRestartFileFormat%isBinary) then
                write(spec%restartFile%unit) meanAccRateSinceStart
            elseif (present(proposal)) then
                ! The extra components are all informational for the end user, otherwise not needed for a restart.
                write(spec%restartFile%unit, spec%restartFile%format) "meanAcceptanceRateSinceStart" & ! LCOV_EXCL_LINE
                                                                    , meanAccRateSinceStart & ! LCOV_EXCL_LINE
                                                                    , "uniqueStateVisitCount" & ! LCOV_EXCL_LINE
                                                                    , proposal%sampleSizeOld & ! LCOV_EXCL_LINE
                                                                    , "proposalAdaptiveScaleSq" & ! LCOV_EXCL_LINE
                                                                    , proposal%scaleSq%running * proposal%scaleSq%default & ! LCOV_EXCL_LINE
                                                                    , "proposalCovLogVol" & ! LCOV_EXCL_LINE
                                                                    , proposal%logSqrtDetOld & ! LCOV_EXCL_LINE
                                                                    , "proposalMean" & ! LCOV_EXCL_LINE
                                                                    , proposal%meanOld & ! LCOV_EXCL_LINE
                                                                    , "proposalChol" & ! LCOV_EXCL_LINE
                                                                    , ((proposal%covLowChoUpp(idim, jdim, 0), idim = 1, jdim - 1), jdim = 2, spec%ndim%vp1)
                                                                    ! Only the upper-triangle of the Cholesky factorization with no delayed rejection.
                                                                    ! We pass Cholesky factorization because the covariance matrix tends to become
                                                                    ! singular in the post-processing step in higher-level languages.
            else
                ! The extra components are all informational for the end user, otherwise not needed for a restart.
                write(spec%restartFile%unit, spec%restartFile%format) "meanAcceptanceRateSinceStart" & ! LCOV_EXCL_LINE
                                                                    , meanAccRateSinceStart
            end if
            flush(spec%restartFile%unit)
        end subroutine writeRestart

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Return a newly sampled state to the kernel routine for acceptance check.
        ! ATTN: The colon index in place of 1 : ndim below avoids the temporary array creation.
        ! The purity breaks only because of the warn() call.
        subroutine setProposalStateNew(spec, proposal, delayedRejectionStage, stateOld, stateNew, rng, err)
            use pm_distMultiNorm, only: setMultiNormRand, uppDia
            use pm_distUnif, only: xoshiro256ssw_type
            use pm_distUnifEll, only: setUnifEllRand
            use pm_ellipsoid, only: isMemberEll
            type(err_type), intent(inout) :: err
            type(spec_type), intent(inout) :: spec
            type(proposal_type), intent(inout) :: proposal
            type(xoshiro256ssw_type), intent(inout) :: rng
            real(RKG), intent(in), contiguous :: stateOld(:)
            real(RKG), intent(out), contiguous :: stateNew(:)
            integer(IK), intent(in) :: delayedRejectionStage
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@setProposalStateNew()"
            character(*, SK), parameter :: WARNING_MSG = PROCEDURE_NAME//SKG_"Number of proposed states out of the objective function domain without any acceptance: "
            character(*, SK), parameter :: FATAL_PREFIX = PROCEDURE_NAME//SKG_"Number of proposed states out of the objective function domain without any acceptance has reach the maximum specified value (domainErrCountMax = "
            character(*, SK), parameter :: FATAL_SUFFIX = SKG_"). This can occur if the domain of the density function is incorrectly specified or if the density function definition implementation or definition is flawed."
            integer(IK) :: domainCheckCounter
            domainCheckCounter = 1_IK
            ! Check for the support Region consistency:
            loopBoundaryCheck: do
                ! \todo
                ! There is an extremely low but non-zero likelihood that the resulting `stateNew` would be the same as `stateOld`.
                ! This could later become problematic because the same state is accepted as new state, potentially corrupting
                ! the proposal covariance computation based on the newly-collected samples.
                if (spec%proposal%is%normal) then
                    call setMultiNormRand(rng, stateNew, stateOld, proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, delayedRejectionStage), uppDia)
                elseif (spec%proposal%is%uniform) then
                    call setUnifEllRand(rng, stateNew, stateOld, proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, delayedRejectionStage), uppDia)
                else
                    error stop "Internal error occurred: proposal distribution unrecognized."
                end if
                ! check if proposal is in the domain.
                if (spec%domain%isFinite) then
                    if (spec%domain%isCube) then
                        if (all(spec%domainCubeLimitLower%val <= stateNew .and. stateNew <= spec%domainCubeLimitUpper%val)) return
                    elseif (spec%domain%isBall) then
                        if (isMemberEll(spec%domainBallCov%inv, spec%domainBallAvg%val, stateNew)) return
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
        !>  \param[in]      sampleState                     :   The two dimensional `(ndim, nsam)` array of accepted states.
        !>  \param[in]      sampleWeight                    :   The one dimensional integer array of the weights of the sampled states in the sampleState.
        !>  \param[in]      proposalAdaptationGreedyEnabled :   The logical flag indicating whether the adaptation is in greedy mode..
        !>  \param[inout]   meanAccRateSinceStart           :   The mean acceptance rate since the start of the sampling.
        !>  \param[out]     proposalAdaptationSampleUsed     :   The output logical flag indicating whether the sampling succeeded.
        !>  \param[out]     proposalAdaptation              :   The output real number in the range `[0,1]` indicating the amount of adaptation,
        !>                                                      with zero indicating no adaptation and one indicating extreme adaptation to the extent
        !>                                                      that the new adapted proposal distribution is completely different from the previous proposal.
        !>  \warning
        !>  This routine must be exclusively called by the leader images.
        !>
        !>  \note
        !>  Other than the call to `warn%show()`, this procedure is `pure`.<br>
        subroutine setProposalAdapted   ( spec & ! LCOV_EXCL_LINE
                                        , proposal & ! LCOV_EXCL_LINE
                                        , sampleState & ! LCOV_EXCL_LINE
                                        , sampleWeight & ! LCOV_EXCL_LINE
                                        , proposalAdaptationGreedyEnabled & ! LCOV_EXCL_LINE
                                        , meanAccRateSinceStart & ! LCOV_EXCL_LINE
                                        , proposalAdaptationSampleUsed & ! LCOV_EXCL_LINE
                                        , proposalAdaptation & ! LCOV_EXCL_LINE
                                        , err & ! LCOV_EXCL_LINE
                                        )

            use pm_sampleMean, only: setMean
            use pm_sampleMean, only: setMeanMerged
            use pm_sampleCov, only: setCovMeanUpdated
            use pm_matrixTrace, only: getMatMulTraceLog
            use pm_sampleCov, only: setCov, setCovMeanMerged
            use pm_matrixChol, only: setMatChol, lowDia, transHerm
            character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@setProposalAdapted()"
            type(spec_type), intent(inout) :: spec
            type(proposal_type), intent(inout) :: proposal
            real(RKG), intent(in), contiguous :: sampleState(:,:)
            integer(IK), intent(in), contiguous :: sampleWeight(:)
            logical(LK), intent(in) :: proposalAdaptationGreedyEnabled
            real(RKG), intent(inout) :: meanAccRateSinceStart ! is intent(out) in restart mode, intent(in) in fresh mode.
            logical(LK), intent(out) :: proposalAdaptationSampleUsed
            real(RKG), intent(out) :: proposalAdaptation
            type(err_type), intent(inout) :: err
            real(RKG) :: sampleStateMeanNew(size(sampleState, 1, IK))
            real(RKG) :: sampleStateMeanCurrent(size(sampleState, 1, IK))
            real(RKG) :: sampleStateCovCholOld(size(sampleState, 1, IK), size(sampleState, 1, IK) + 1)
            real(RKG) :: sampleStateCovLowCurrent(size(sampleState, 1, IK), size(sampleState, 1, IK) + 1)
            real(RKG) :: logSqrtDetNew, logSqrtDetSum, adaptiveScale, adaptiveScaleSq, fracA
            integer(IK) :: sampleSizeOld, sampleSizeCurrent
            integer(IK):: jdim, nsam, info
            integer(IK), parameter :: dim = 2_IK
            logical(LK) :: proposalScalingNeeded
            logical(LK) :: nonSingularSample

            nsam = size(sampleState, 2, IK)
            sampleSizeOld = proposal%sampleSizeOld ! this is kept only for restoration of proposal%sampleSizeOld, if needed.
            proposalScalingNeeded = spec%targetAcceptanceRate%enabled

            ! First if there are less than spec%ndim%vp1 points for new covariance computation, then just scale the covariance and return.

            ! This will be used to recover the old covariance in case of proposal update failure, and to compute the adaptation measure.
            sampleStateCovCholOld = proposal%covLowChoUpp(:, :, 0)

            nonSingularSample = spec%ndim%val < nsam
           !proposalAdaptationSampleUsed = nonSingularSample ! update proposal only if sample is non-singular.
            proposalAdaptationSampleUsed = nonSingularSample .or. 1_IK < proposal%sampleSizeOld ! update sample for any input singular or non-singular sample.

            mergeCovMat_block: if (proposalAdaptationSampleUsed) then

                ! Get the upper covariance matrix and mean of the new sample.

                if (proposalAdaptationGreedyEnabled) then
                    sampleSizeCurrent = nsam
                    call setMean(sampleStateMeanCurrent, sampleState, dim)
                    if (nonSingularSample) call setCov(sampleStateCovLowCurrent(:, 1 : spec%ndim%val), lowDia, sampleStateMeanCurrent, sampleState, dim)
                else
                    call setMean(sampleStateMeanCurrent, sampleState, dim, sampleWeight, sampleSizeCurrent)
                    if (nonSingularSample) call setCov(sampleStateCovLowCurrent(:, 1 : spec%ndim%val), lowDia, sampleStateMeanCurrent, sampleState, dim, sampleWeight, sampleSizeCurrent)
                end if

                ! Combine old and new covariance matrices if both exist.

                if (1_IK < proposal%sampleSizeOld) then

                    fracA = real(proposal%sampleSizeOld, RKG) / real(proposal%sampleSizeOld + sampleSizeCurrent, RKG)

                    nonSingular_block: if (nonSingularSample) then

                        ! First scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor.

                        if (proposal%scaleSq%default /= 1._RKG) then
                            ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                            !do concurrent(jdim = 1 : spec%ndim%val)
                            do jdim = 1, spec%ndim%val
                                sampleStateCovLowCurrent(jdim : spec%ndim%val, jdim) = sampleStateCovLowCurrent(jdim : spec%ndim%val, jdim) * proposal%scaleSq%default
                            end do
                        end if

                        ! Now combine it with the old covariance matrix.
                        ! Do not set the full boundaries' range `(1 : spec%ndim%val)` for `proposal%covLowChoUpp` in the following subroutine call.
                        ! Setting the boundaries forces the compiler to generate a temporary array.

                        call setCovMeanMerged   ( proposal%covLowChoUpp(:, 1 : spec%ndim%val, 0) & ! LCOV_EXCL_LINE
                                                , sampleStateMeanNew & ! LCOV_EXCL_LINE
                                                , sampleStateCovLowCurrent(:, 1 : spec%ndim%val) & ! LCOV_EXCL_LINE
                                                , sampleStateMeanCurrent & ! LCOV_EXCL_LINE
                                                , sampleStateCovCholOld(:, 1 : spec%ndim%val) & ! LCOV_EXCL_LINE
                                                , proposal%meanOld & ! LCOV_EXCL_LINE
                                                , fracA & ! LCOV_EXCL_LINE
                                                , lowDia & ! LCOV_EXCL_LINE
                                                )
                        proposal%meanOld(1 : spec%ndim%val) = sampleStateMeanNew

                    else nonSingular_block ! singular new sample.

                        call setCovMeanUpdated  ( covA = proposal%covLowChoUpp(:, 1 : spec%ndim%val, 0) & ! LCOV_EXCL_LINE
                                                , meanA = proposal%meanOld & ! LCOV_EXCL_LINE
                                                , meanB = sampleStateMeanCurrent & ! LCOV_EXCL_LINE
                                                , fracA = fracA & ! LCOV_EXCL_LINE
                                                , subset = lowDia & ! LCOV_EXCL_LINE
                                                )

                    end if nonSingular_block

                else ! first occurrence of nonSingularSample.

                    ! There is no prior old Covariance matrix to combine with the new one from the new sampleState

                    proposal%meanOld(1 : spec%ndim%val) = sampleStateMeanCurrent
                    proposal%sampleSizeOld = sampleSizeCurrent

                    ! Copy and then scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor.

                    if (proposal%scaleSq%default == 1._RKG) then
                        ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                        !do concurrent(jdim = 1 : spec%ndim%val)
                        do jdim = 1, spec%ndim%val
                            proposal%covLowChoUpp(jdim : spec%ndim%val, jdim, 0) = sampleStateCovLowCurrent(jdim : spec%ndim%val, jdim)
                        end do
                    else
                        do jdim = 1, spec%ndim%val
                            proposal%covLowChoUpp(jdim : spec%ndim%val, jdim, 0) = sampleStateCovLowCurrent(jdim : spec%ndim%val, jdim) * proposal%scaleSq%default
                        end do
                    end if

                end if

                ! Now get the Cholesky factorization.
                ! WARNING: Do not set the full boundaries' range `(1 : spec%ndim%val)` for the first index of `proposal%covLowChoUpp` in the following subroutine call.
                ! WARNING: Setting the boundaries forces the compiler to generate a temporary array.
                ! The `:` instead of `1 : spec%ndim%val` avoids temporary array creation.

                call setMatChol(proposal%covLowChoUpp(:, 1 : spec%ndim%val, 0), lowDia, info, proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, 0), transHerm)

                proposalAdaptationSampleUsed = info == 0_IK
                posDefCheck_block: if (proposalAdaptationSampleUsed) then

                    proposal%sampleSizeOld = proposal%sampleSizeOld + sampleSizeCurrent

                else posDefCheck_block ! singularityOccurred

                    proposalAdaptation = 0._RKG
                    proposal%sampleSizeOld = sampleSizeOld

                    ! It may be a good idea to add a warning message printed out here for the singularity occurrence.
                    call spec%disp%warn%show(SKG_"Singularity occurred while adapting the proposal distribution covariance matrix.", tmsize = 0_IK, bmsize = 0_IK)

                    ! Recover the old covariance matrix + Cholesky factor.
                    proposal%covLowChoUpp(:, :, 0) = sampleStateCovCholOld

                end if posDefCheck_block

            else mergeCovMat_block ! Happens only if this is the first update and the first sample is singular.

                ! If the first covariance merging has not occurred yet,
                ! set the scaling factor appropriately to shrink the covariance matrix.

                !proposalScalingNeeded = proposalScalingNeeded .or. 1_IK == proposal%sampleSizeOld ! scale only when sample is singular in the first try.
                proposalScalingNeeded = .true._LK

            end if mergeCovMat_block

            !! perform global adaptive scaling is requested
            !if (spec%targetAcceptanceRate%enabled) then
            !    if (meanAccRateSinceStart < spec%targetAcceptanceRate%aim) then
            !        proposal%scaleSq%running = mc_maxScaleFactorSq**((spec%targetAcceptanceRate%aim-meanAccRateSinceStart)/spec%targetAcceptanceRate%aim)
            !    else
            !        proposal%scaleSq%running = mc_maxScaleFactorSq**((spec%targetAcceptanceRate%aim-meanAccRateSinceStart)/(1._RKG-spec%targetAcceptanceRate%aim))
            !    end if
            !    proposal%scaleSq%running = proposal%scaleSq%running * (meanAccRateSinceStart/spec%targetAcceptanceRate%aim)**spec%ndim%inv
            !end if

            ! Read the restart file (in dry mode).

            if (spec%image%is%leader .and. spec%run%is%dry) then
                if (proposalAdaptationSampleUsed .or. proposalScalingNeeded) then
                    call readRestart(spec, meanAccRateSinceStart)
                else
                    call readRestart(spec, meanAccRateSinceStart, nskip = 0_IK)
                end if
            end if

            ! Adjust the scale of the covariance matrix and the Cholesky factor, if needed.

            proposalScaling_block: if (proposalScalingNeeded) then

                if (spec%targetAcceptanceRate%enabled) then
                    if (meanAccRateSinceStart < spec%targetAcceptanceRate%val(1) .or. spec%targetAcceptanceRate%val(2) < meanAccRateSinceStart) then
                        adaptiveScale = (meanAccRateSinceStart / spec%targetAcceptanceRate%aim) ** spec%ndim%invhalf
                        proposal%scaleSq%running = adaptiveScale ** 2
                        adaptiveScaleSq = proposal%scaleSq%running
                        !block
                        !    use pm_distUnif, only: getUnifRand
                        !    integer, save :: counter = 0_IK
                        !    counter = counter - 1
                        !    proposal%scaleSq%running = proposal%scaleSq%running * exp(-counter*getUnifRand(-1._RKG, 1._RKG)/1.e4_RK)
                        !    !use pm_statistics, only: getUnifRand
                        !    !proposal%scaleSq%running = proposal%scaleSq%running * exp(real(getUnifRand(-1_IK, 1_IK), RK))
                        !    !write(*,*) counter, proposal%scaleSq%running
                        !end block
                    end if
                else
                    ! This should happen only for the first update with singular samples.
                    adaptiveScale = 0.5_RKG
                    adaptiveScaleSq = 0.25_RKG
                end if

                proposal%covLowChoUpp(1 : spec%ndim%val, 1, 0) = proposal%covLowChoUpp(1 : spec%ndim%val, 1, 0) * adaptiveScaleSq ! Update first column of covariance matrix.
                ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                !do concurrent(jdim = 2 : spec%ndim%val)
                do jdim = 2, spec%ndim%val
                    proposal%covLowChoUpp(1 : jdim - 1, jdim, 0) = proposal%covLowChoUpp(1 : jdim - 1, jdim, 0) * adaptiveScale ! Update the Cholesky factor.
                    proposal%covLowChoUpp(jdim : spec%ndim%val, jdim, 0) = proposal%covLowChoUpp(jdim : spec%ndim%val, jdim, 0) * adaptiveScaleSq ! Update covariance matrix.
                end do
                proposal%covLowChoUpp(1 : spec%ndim%val, spec%ndim%vp1, 0) = proposal%covLowChoUpp(1 : spec%ndim%val, spec%ndim%vp1, 0) * adaptiveScale ! Update last column of Cholesky factor.

            end if proposalScaling_block

            ! Compute the adaptivity only if any updates has occurred.

            proposalAdaptationEstimate_block: if (proposalAdaptationSampleUsed .or. proposalScalingNeeded) then

                logSqrtDetNew = getMatMulTraceLog(proposal%covLowChoUpp(:, 2 : spec%ndim%vp1, 0))
                ! Use the universal upper bound.
                !block
                !use pm_io
                !call disp%show(sampleStateCovLowCurrent)
                !end block
                ! The `do concurrent` version is problematic for OpenMP builds of the ParaMonte library with the Intel ifort 2021.8 on heap, yielding segmentation fault.
                !do concurrent(jdim = 1 : spec%ndim%val)
                do jdim = 1, spec%ndim%val
                    sampleStateCovLowCurrent(jdim : spec%ndim%val, jdim) = 0.5_RKG * (proposal%covLowChoUpp(jdim : spec%ndim%val, jdim, 0) + sampleStateCovCholOld(jdim : spec%ndim%val, jdim)) ! dummy
                end do
                call setMatChol(sampleStateCovLowCurrent(:, 1 : spec%ndim%val), lowDia, info, sampleStateCovLowCurrent(:, 2 : spec%ndim%vp1), transHerm)
                if (info == 0_IK) then
                    logSqrtDetSum = getMatMulTraceLog(sampleStateCovLowCurrent(:, 2 : spec%ndim%vp1))
                else ! singularityOccurred ! LCOV_EXCL_START
                    err%occurred = .true._LK
                    err%msg = NL1//PROCEDURE_NAME//SKG_"@line::"//getStr(__LINE__)//SKG_": "//&
                    SKG_"Singular lower covariance matrix detected at column ("//getStr(info)//SKG_") while computing the proposalAdaptation measure:"//NL2
                    do info = 1, size(sampleStateCovLowCurrent, 1, IK)
                        err%msg = err%msg//getStr(sampleStateCovLowCurrent(info,:))//NL1
                    end do
                    err%msg = err%msg//NL1//&
                    SKG_": Error occurred while computing the Cholesky factorization of a matrix needed for the computation of the adaptation measure. &
                    &Such error is highly unusual, and requires an in depth investigation of the case. &
                    &It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations. &
                    &In such case, consider increasing the value of the input variable proposalAdaptationPeriod. &
                    &It may also be that your input objective function has been incorrectly implemented. &
                    &Also, ensure a correct value of `ndim` to the ParaMonte sampler routine, representing the domain size. &
                    &Otherwise, restarting the simulation might resolve the error. If the error persists, please inform the developers at: "//NL2//PARAMONTE_WEB_ISSUES
                    return
                end if ! LCOV_EXCL_STOP
                proposalAdaptation = 1._RKG - exp(0.5_RKG * (proposal%logSqrtDetOld + logSqrtDetNew) - logSqrtDetSum) ! totalVariationUpperBound
                if (0._RKG < proposalAdaptation) then
                    proposalAdaptation = sqrt(proposalAdaptation) ! totalVariationUpperBound
                elseif (proposalAdaptation < 0._RKG) then
                    spec%msg = PROCEDURE_NAME//SKG_": Non-positive adaptation measure detected" ! LCOV_EXCL_LINE
                    if (abs(proposalAdaptation) < sqrt(epsilon(0._RKG))) spec%msg = spec%msg//SKG_" (possibly due to round-off error)" ! LCOV_EXCL_LINE
                    call spec%disp%warn%show(spec%msg//SKG_": "//getStr(proposalAdaptation)) ! LCOV_EXCL_LINE
                    proposalAdaptation = 0._RKG ! LCOV_EXCL_LINE
                end if
                proposal%logSqrtDetOld = logSqrtDetNew
                !block
                !integer, save :: counter = 0
                !counter = counter + 1
                !!if (counter==1) then
                !if (1._RKG < proposalAdaptation) then
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

            else proposalAdaptationEstimate_block

                proposalAdaptation = 0._RKG

            end if proposalAdaptationEstimate_block

            ! Read the restart file (in new mode).

            if (spec%image%is%leader .and. spec%run%is%new) then
                ! This extra check avoids redundant restart writes when no scaling has occurred.
                if (proposalAdaptationSampleUsed .or. proposalScalingNeeded) then
                    call writeRestart(spec, meanAccRateSinceStart, proposal)
                else
                    call writeRestart(spec, meanAccRateSinceStart)
                end if
            end if
            !if (spec%image%is%leader .and. spec%run%is%new) call writeRestart(spec, proposal, meanAccRateSinceStart)

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
