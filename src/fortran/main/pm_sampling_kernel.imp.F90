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

#if     TESTING_ENABLED && (CAF_ENABLED || MPI_ENABLED)
#define FAILED_IMAGE(FAILED)isFailedImage(FAILED)
#else
#define FAILED_IMAGE(FAILED)FAILED
#endif
#if     ParaDISE_ENABLED
#define SET_ParaDISE(X)X
#else
#define SET_ParaDISE(X)
#endif
#if     CAF_ENABLED || MPI_ENABLED
#define SET_CAFMPI(X)X
#define SET_SERIAL(X)
#else
#define SET_CAFMPI(X)
#define SET_SERIAL(X)X
#endif
#if     CAF_ENABLED
#define SET_CAF(X)X
#else
#define SET_CAF(X)
#endif
#if     MPI_ENABLED
        use mpi !mpi_f08, only: mpi_bcast, mpi_gather, mpi_allreduce, mpi_sum, mpi_integer, mpi_double_precision, mpi_comm_world
#define SET_MPI(X)X
#else
#define SET_MPI(X)
#endif
#if     OMP_ENABLED
#define SET_OMP(X)X
#else
#define SET_OMP(X)
#endif
#if     DEBUG_ENABLED
#define SET_DEBUG(X)X
#else
#define SET_DEBUG(X)
#endif
        ! \bug
        ! Avoid Intel ifort bug for too many `use` statements in a submodule by placing them all in the submodule header.
#if     CHECK_ENABLED
        use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG)call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
        ! Define error handling.
#define RETURN_IF_FAILED(LINE,FAILED,MSG) \
if (FAILED) then; \
err%occurred = .true._LK; \
err%msg = PROCEDURE_NAME//getFine(__FILE__, LINE)//SKG_": "//trim(MSG); \
call spec%disp%stop%show(err%msg); \
return; \
end if;
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     ParaDISE_ENABLED || ParaDRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_err, only: getFine
        use pm_err, only: err_type
        use pm_val2str, only: getStr
        use pm_timer, only: timer_type
        use pm_kind, only: SKG => SK, SK, IK, LK
        use pm_paramonte, only: PARAMONTE_WEB_ISSUES
        use pm_parallelism, only: isFailedImage
        use pm_sampling_proposal, only: NL1, NL2
        use pm_sampling_proposal, only: spec_type
        use pm_sampling_proposal, only: stat_type
        use pm_sampling_proposal, only: getErrChainRead
        use pm_sampling_proposal, only: getErrChainWrite
        use pm_sampling_proposal, only: chainFileColName
        use pm_sampling_proposal, only: isFailedChainResize
        use pm_sampling_proposal, only: proposal_type, setProposalAdapted, readRestart, writeRestart, setProposalStateNew
#if     CAF_ENABLED || MPI_ENABLED
        use pm_sampling_proposal, only: bcastProposalAdaptation
#endif
        ! used by the kernel routines.
        use pm_distUnif, only: setUnifRand
        use pm_arrayResize, only: setResized
        use pm_sysPath, only: isFailedRemove
        use pm_mathLogSubExp, only: getLogSubExp
        use pm_io, only: setContentsFrom, setContentsTo
        use, intrinsic :: iso_c_binding, only: c_sizeof
        use pm_kind, only: RKD

        implicit none
#if     ParaDISE_ENABLED
        character(*,SKG), parameter :: MODULE_NAME = SK_"@pm_sampling_kernel_dise"
#elif   ParaDRAM_ENABLED
        character(*,SKG), parameter :: MODULE_NAME = SK_"@pm_sampling_kernel_dram"
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        abstract interface
#if         OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
            function getLogFunc_proc(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
                import :: RKG, RKD
                real(RKG), intent(inout), contiguous :: logFuncState(0:, :)
                real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
                real(RKG) :: mold
            end function
#else
            function getLogFunc_proc(state) result(logFunc)
                import :: RKG
                real(RKG), intent(in), contiguous :: state(:)
                real(RKG) :: logFunc
            end function
#endif
        end interface

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ParaDISE_ENABLED
        subroutine killMeAlreadyCMake_RK5(); use pm_sampling_proposal_dise_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK4(); use pm_sampling_proposal_dise_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK3(); use pm_sampling_proposal_dise_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK2(); use pm_sampling_proposal_dise_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK1(); use pm_sampling_proposal_dise_RK1, only: RKG; end subroutine
#elif   ParaDRAM_ENABLED
        subroutine killMeAlreadyCMake_RK5(); use pm_sampling_proposal_dram_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK4(); use pm_sampling_proposal_dram_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK3(); use pm_sampling_proposal_dram_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK2(); use pm_sampling_proposal_dram_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK1(); use pm_sampling_proposal_dram_RK1, only: RKG; end subroutine
#elif   ParaNest_ENABLED
        subroutine killMeAlreadyCMake_RK5(); use pm_sampling_proposal_nest_RK5, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK4(); use pm_sampling_proposal_nest_RK4, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK3(); use pm_sampling_proposal_nest_RK3, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK2(); use pm_sampling_proposal_nest_RK2, only: RKG; end subroutine
        subroutine killMeAlreadyCMake_RK1(); use pm_sampling_proposal_nest_RK1, only: RKG; end subroutine
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !>  Return an error if the sampler kernel run failed.
        !>
        !>  \warning
        !>  This procedure changes the components of the global-scope variables `spec` and `proposal`. Otherwise, it is pure.
        function getErrKernelRun(getLogFunc, spec, stat, proposal) result(err)
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: getErrKernelRun
#endif
            procedure(getLogFunc_proc)          :: getLogFunc
            type(proposal_type) , intent(inout) :: proposal
            type(spec_type)     , intent(inout) :: spec
            type(stat_type)     , intent(inout) :: stat
            type(err_type)                      :: err
#include    "pm_sampling_kernel.inc.F90"
        end function getErrKernelRun

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !>  Write the latest simulated state to the output chain file contents.
        subroutine writeCFC(spec, stat, proposalAdaptation)
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: writeCFC
#endif
            type(spec_type), intent(in) :: spec
            type(stat_type), intent(in) :: stat
            real(RKG), intent(in), contiguous :: proposalAdaptation(:)
            integer(IK) :: j
            ! if new point has been sampled, write the previous sampled point to output file.
            blockOutputWrite: if (0_IK < stat%numFunCallAccepted) then
                if (spec%outputChainFileFormat%isCompact) then
                    write(spec%chainFile%unit,spec%chainFile%format%rows) stat%cfc%processID            (stat%numFunCallAccepted) &
                                                                        , stat%cfc%delayedRejectionStage(stat%numFunCallAccepted) &
                                                                        , stat%cfc%meanAcceptanceRate   (stat%numFunCallAccepted) &
                                                                        , stat%cfc%proposalAdaptation   (stat%numFunCallAccepted) &
                                                                        , stat%cfc%burninLocation       (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleWeight         (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleLogFunc        (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleState          (1 : spec%ndim%val, stat%numFunCallAccepted)
                elseif (spec%outputChainFileFormat%isBinary) then
                    write(spec%chainFile%unit                           ) stat%cfc%processID            (stat%numFunCallAccepted) &
                                                                        , stat%cfc%delayedRejectionStage(stat%numFunCallAccepted) &
                                                                        , stat%cfc%meanAcceptanceRate   (stat%numFunCallAccepted) &
                                                                        , stat%cfc%proposalAdaptation   (stat%numFunCallAccepted) &
                                                                        , stat%cfc%burninLocation       (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleWeight         (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleLogFunc        (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleState          (1 : spec%ndim%val, stat%numFunCallAccepted)
                elseif (spec%outputChainFileFormat%isVerbose) then
                    do j = 1, stat%cfc%sampleWeight(stat%numFunCallAccepted)
                    write(spec%chainFile%unit,spec%chainFile%format%rows) stat%cfc%processID            (stat%numFunCallAccepted) &
                                                                        , stat%cfc%delayedRejectionStage(stat%numFunCallAccepted) &
                                                                        , stat%cfc%meanAcceptanceRate   (stat%numFunCallAccepted) &
                                                                        , proposalAdaptation(j) &
                                                                       !, stat%cfc%proposalAdaptation(stat%numFunCallAccepted) &
                                                                        , stat%cfc%burninLocation       (stat%numFunCallAccepted) &
                                                                        , 1_IK &
                                                                        , stat%cfc%sampleLogFunc        (stat%numFunCallAccepted) &
                                                                        , stat%cfc%sampleState          (1 : spec%ndim%val, stat%numFunCallAccepted)
                    end do
                end if
            end if blockOutputWrite
        end subroutine writeCFC

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !> \brief
        !> Update the user with the simulation progress information.
        !>
        !> \remark
        !> Objects that change state in this subroutine are: `timer`, `stat%progress%timeElapsedSinceStartInSeconds`, `stat%numFunCallAcceptedLastReport`.
        subroutine reportProgress(spec, stat, timeLeft)!, sumAccrAccRej
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: reportProgress
#endif
            use pm_strASCII, only: CR
            use pm_arraySplit, only: setSplit
            use iso_fortran_env, only: output_unit
            class(spec_type), intent(inout) :: spec
            class(stat_type), intent(inout) :: stat
            real(RKD), intent(in), optional :: timeLeft
            real(RKD) :: timeElapsedSinceLastReportInSeconds
            real(RKD) :: timeRemainedToFinishInSeconds
            real(RKG) :: meanAccRateSinceLastReport
            real(RKG) :: meanAccRateSinceStart
            character(1 + len(spec%reportFile%indent) + 3 * (25 + 3), SK) :: status
            if (spec%run%is%new) then
                stat%timer%clock = stat%timer%time()
                if (stat%progress%clock < 0._RKG) stat%progress%clock = stat%timer%start
                timeElapsedSinceLastReportInSeconds = stat%timer%clock - stat%progress%clock
                stat%progress%timeElapsedSinceStartInSeconds = stat%progress%timeElapsedSinceStartInSeconds + timeElapsedSinceLastReportInSeconds
                stat%progress%clock = stat%timer%clock
                if (present(timeLeft)) then
                    timeRemainedToFinishInSeconds = timeLeft
                else
                    timeRemainedToFinishInSeconds = stat%progress%timeElapsedSinceStartInSeconds * real(spec%outputChainSize%val - stat%numFunCallAccepted, RKD) / stat%numFunCallAccepted
                end if
                CHECK_ASSERTION(__LINE__, stat%progress%counterPRP == spec%outputReportPeriod%val .or. spec%outputReportPeriod%val == 0_IK .or. timeRemainedToFinishInSeconds == 0._RKG, \
                SK_"The condition `counterPRP == outputReportPeriod .or. outputReportPeriod == 0 .or. timeLeft == 0` must hold. counterPRP, outputReportPeriod, timeLeft = "//\
                getStr([stat%progress%counterPRP, spec%outputReportPeriod%val])//SK_", "//getStr(timeRemainedToFinishInSeconds))
                meanAccRateSinceStart = real(stat%numFunCallAccepted, RKG) / real(stat%numFunCallAcceptedRejected, RKG)
                meanAccRateSinceLastReport = (stat%numFunCallAccepted - stat%numFunCallAcceptedLastReport) * spec%outputReportPeriod%inv
                stat%numFunCallAcceptedLastReport = stat%numFunCallAccepted
                write(spec%progressFile%unit, spec%progressFile%format%rows ) stat%numFunCallAcceptedRejected &
                                                                            , stat%numFunCallAccepted &
                                                                            , meanAccRateSinceStart &
                                                                            , meanAccRateSinceLastReport &
                                                                            , timeElapsedSinceLastReportInSeconds &
                                                                            , stat%progress%timeElapsedSinceStartInSeconds &
                                                                            , timeRemainedToFinishInSeconds
                flush(spec%progressFile%unit)
            else
                block
                    integer(IK) :: lenrec
                    character(500, SK) :: record
                    integer(IK), allocatable :: index(:,:)
                    integer(IK) :: numFunCallAcceptedRejectedLastReport
                    read(spec%progressFile%unit, "(A)") record
                    lenrec = len_trim(record, IK)
                    call setSplit(index, record(1:lenrec), spec%outputSeparator%val)
                    read(record(index(1,1) : index(2,1)), *) numFunCallAcceptedRejectedLastReport
                    read(record(index(1,2) : index(2,2)), *) stat%numFunCallAcceptedLastReport
                    read(record(index(1,3) : index(2,3)), *) meanAccRateSinceStart
                    read(record(index(1,4) : index(2,4)), *) meanAccRateSinceLastReport
                    read(record(index(1,5) : index(2,5)), *) timeElapsedSinceLastReportInSeconds
                    read(record(index(1,6) : index(2,6)), *) stat%progress%timeElapsedSinceStartInSeconds
                    read(record(index(1,7) : index(2,7)), *) timeRemainedToFinishInSeconds ! estimatedTimeToFinishInSeconds
                end block
                stat%progress%clock = stat%timer%time()
            end if
            ! report progress in the standard output
            if (spec%image%is%first) then
                write(status, fmt = "(A1,A,*(A25,3X))") CR, spec%reportFile%indent, & ! LCOV_EXCL_LINE
                getStr(getStr(stat%numFunCallAccepted)//SK_" / "//getStr(stat%numFunCallAcceptedRejected), SK_"(125A)"), & ! LCOV_EXCL_LINE
                getStr(trim(adjustl(getStr(meanAccRateSinceLastReport, SK_"(1F11.4)")))//SK_" / "//trim(adjustl(getStr(stat%numFunCallAccepted / real(stat%numFunCallAcceptedRejected), SK_"(1F11.4)"))), SK_"(125A)"), & ! LCOV_EXCL_LINE
                getStr(trim(adjustl(getStr(stat%progress%timeElapsedSinceStartInSeconds, SK_"(1F11.4)")))//SK_" / "//trim(adjustl(getStr(min(999999._RKD, timeRemainedToFinishInSeconds), SK_"(1F11.4)"))), SK_"(125A)")
                call spec%disp%show(status, unit = output_unit, advance = SK_"NO", tmsize = 0_IK, bmsize = 0_IK)
                if (timeRemainedToFinishInSeconds == 0._RKG) call spec%disp%skip(unit = output_unit)
                !call execute_command_line(" ")
                flush(output_unit)
                ! LCOV_EXCL_STOP
            end if
            !numFunCallAcceptedRejectedLastReport = stat%numFunCallAcceptedRejected
            stat%progress%counterPRP = 0_IK
        end subroutine reportProgress

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine initProgressReport(spec)
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: initProgressReport
#endif
            use iso_fortran_env, only: output_unit
            type(spec_type), intent(inout) :: spec
            character(:, SK), allocatable :: header
            integer(IK), parameter :: SEGLEN = 25
            if (spec%image%is%first) then
                ! LCOV_EXCL_START
                header = spec%reportFile%indent//"Accepted/Total Func. Call   "//"Dynamic/Overall Acc. Rate   "//"Elapsed/Remained Time [s] "//NL1//&
                repeat(" ",len(spec%reportFile%indent))//repeat("=",SEGLEN)//"   "//repeat("=",SEGLEN)//"   "//repeat("=",SEGLEN)//" "
                call spec%disp%show(header, unit = output_unit, format = SK_"(A)", tmsize = 0_IK, bmsize = 0_IK)
                !call execute_command_line(" ")
                deallocate(header)
                flush(output_unit)
                ! LCOV_EXCL_STOP
            end if
        end subroutine initProgressReport

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !> Return the predicted burnin location within the currently sampled chain.
        !>
        !> \param[in]   lenLogFunc  :   The length of the input `lenLogFunc`.
        !> \param[in]   refLogFunc  :   The reference logFunc value with respect to which the relative importance of the points is gauged.
        !> \param[in]   logFunc     :   The vector of length `lenLogFunc` of log(function) values.
        !>
        !>  \devnote
        !>  This is a `simple` function.
        !>
        !> \return
        !> `burninLoc` : The location of burnin in the input `logFunc` vector.
        pure function getBurninLoc(lenLogFunc, refLogFunc, logFunc) result(burninLoc)
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: getBurninLoc
#endif
            integer(IK), intent(in) :: lenLogFunc
            real(RKG), intent(in) :: refLogFunc, logFunc(lenLogFunc)
            real(RKG) :: negLogIncidenceProb
            integer(IK) :: burninLoc
            negLogIncidenceProb = log(real(lenLogFunc, RKG))
            burninLoc = 0_IK
            do
                burninLoc = burninLoc + 1_IK
                if (burninLoc < lenLogFunc .and. refLogFunc - logFunc(burninLoc) > negLogIncidenceProb) cycle
                !burninLoc = burninLoc - 1
                exit
            end do
        end function getBurninLoc

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  RETURN_IF_FAILED
#undef  SET_ParaDISE

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  RETURN_IF_FAILED
#undef  SET_CAFMPI
#undef  SET_SERIAL
#undef  SET_DEBUG
#undef  SET_CAF
#undef  SET_MPI
#undef  SET_OMP