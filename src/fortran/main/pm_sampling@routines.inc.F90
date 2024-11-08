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

!>  \brief
!>  This file contains the implementations of the generic interfaces in [pm_sampling]@(ref pm_sampling).
!>
!>  \test
!>  [test_pm_sampling](@ref test_pm_sampling)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RETURN_IF_FAILED(LINE,FAILED,MSG) \
if (FAILED) then; \
err%occurred = .true._LK; \
err%msg = PROCEDURE_NAME//getLine(LINE)//SKG_": "//NL1//MSG; \
call spec%disp%stop%show(err%msg); \
return; \
end if;
        ! Define intel-specific shared property of files on Windows.
#if     INTEL_ENABLED && WINDOWS_ENABLED
#define SHARED, shared
#else
#define SHARED
#endif
        ! Define runtime error handling mode.
#if     TESTING_ENABLED && (CAF_ENABLED || MPI_ENABLED || OMP_ENABLED)
#define FAILED_IMAGE(FAILED)isFailedImage(FAILED)
#else
#define FAILED_IMAGE(FAILED)FAILED
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
!#if     CAF_ENABLED
!#define COINDEX(I)[I]
!#else
!#define COINDEX(I)
!#endif
        ! Set the traceback information.
!#define TRACE(LINE) SKG_"@file::"//getStr(__FILE__)//getLine(LINE)//PROCEDURE_NAME
!        ! Set the alert optional input arguments.
!#define ALERT_OPTIONS unit = reportFileUnit, width = disp%width, prefix = method, tmsize = tmsize, bmsize = bmsize

        !%%%%%%%%%%%%%%%%%%
#if     runParaDRAM_ENABLED
        !%%%%%%%%%%%%%%%%%%

        abstract interface
#if     OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
        recursive function getLogFuncPtr_proc(logFuncState, ndim, njob, avgTimePerFunCall, avgCommPerFunCall) result(mold) bind(C)
            import :: IK, RKG, RKD
            integer(IK), intent(in), value :: ndim, njob
            real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
            real(RKG), intent(inout) :: logFuncState(ndim, njob)
            real(RKG) :: mold
        end function
#else
        recursive function getLogFuncPtr_proc(state, ndim) result(logFunc) bind(C)
            import :: IK, RKG
            integer(IK), intent(in), value :: ndim
            real(RKG), intent(in) :: state(ndim)
            real(RKG) :: logFunc
        end function
#endif
        end interface
        procedure(getLogFuncPtr_proc), pointer :: getLogFuncPtr
        !procedure(real(RKG)), pointer :: getLogFuncPtr
        type(paradram_type) :: sampler
        type(err_type) :: err
        integer(IK) :: lenin
        stat = 0_IK

        ! Reconstruct the input file.

        if (present(input)) then
            lenin = 1
            do
                if (input(lenin) == c_null_char) exit
                lenin = lenin + 1
            end do
            if (1 < lenin) sampler%inputFile = getCharSeq(input(1 : lenin - 1))
        end if

        !if (0 < lenin) then
        !    allocate(character(lenin) :: sampler%inputFile)
        !    do iell = 1, lenin
        !        sampler%inputFile(iell : iell) = input(iell)
        !    end do
        !end if

        ! Associate the input C procedure pointer to a Fortran procedure pointer.

        call c_f_procpointer(cptr = getLogFunc, fptr = getLogFuncPtr)
        err = getErrSampling(sampler, getLogFuncWrapper, ndim) ! run ParaDRAM.
        if (err%occurred) stat = 1_IK
        nullify(getLogFuncPtr)

    contains

#if     OMP_ENABLED && (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED)
        recursive function getLogFuncWrapper(logFuncState, avgTimePerFunCall, avgCommPerFunCall) result(mold)
            real(RKG), intent(inout), contiguous :: logFuncState(0:,:)
            real(RKD), intent(inout) :: avgTimePerFunCall, avgCommPerFunCall
            real(RKG) :: mold
            mold = getLogFuncPtr(logFuncState, ndim, size(logFuncState, 2, IK), avgTimePerFunCall, avgCommPerFunCall)
        end function
#else
        recursive function getLogFuncWrapper(state) result(logFunc)
            real(RKG), intent(in), contiguous :: state(:)
            real(RKG) :: logFunc
            logFunc = getLogFuncPtr(state, ndim)
            !write(*,"(1I5,5E30.20,:,', ')") precision(state), state!, logFunc
        end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getErrSampling_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        use pm_matrixClass, only: isMatClass, posdefmat
        use pm_sampling_kernel, only: getErrChainRead, getErrChainWrite, isFailedChainResize
        use pm_sampling_kernel, only: spec_type, stat_type, proposal_type
        use pm_sampling_kernel, only: chainFileColName
        use pm_sampling_kernel, only: getErrKernelRun
        use pm_sampling_kernel, only: NL1, NL2, RKG
#if     ParaDISE_ENABLED
#define SET_DRAMDISE(X)X
#define SET_ParaNest(X)
        character(*,SKG), parameter :: method = SKG_"ParaDISE"
#elif   ParaDRAM_ENABLED
#define SET_DRAMDISE(X)X
#define SET_ParaNest(X)
        character(*,SKG), parameter :: method = SKG_"ParaDRAM"
#elif   ParaNest_ENABLED
#define SET_ParaNest(X)X
#define SET_DRAMDISE(X)
        character(*,SKG), parameter :: method = SKG_"ParaNest"
#else
#error  "Unrecognized interface."
#endif
        !>  \bug
        !>  \status \unresolved
        !>  \source \ifort{2021.11.1 20231117}
        !>  \desc
        !>  \ifort returns an *already allocated error with the statement `call setResized(scaling, lenScaling)` which persists in both release and debug modes.<br>
        !>  The full debug message is the following:<br>
        !>  \code{.sh}
        !>      forrtl: severe (151): allocatable array is already allocated
        !>      Image              PC                Routine            Line        Source
        !>      libparamonte.so    00007FCFBA641EB8  Unknown               Unknown  Unknown
        !>      libparamonte.so    00007FCFB6E6D4A3  pm_arrayresize_MP         170  pm_arrayResize@routines.inc.F90
        !>      libparamonte.so    00007FCFB72B5BD6  pm_parallelism_MP          77  pm_parallelism@routines.inc.F90
        !>      libparamonte.so    00007FCFB6263F43  pm_sampling_MP_ge         394  pm_sampling@routines.inc.F90
        !>      libparamonte.so    00007FCFB61F0194  runParaDRAMD              136  pm_sampling@routines.inc.F90
        !>  \endcode
        !>  Note that the line numbers for this file in the message above have changed because of code change in this file.<br>
        !>  This error does not occur when the library is compiled with \gfortran{13}.<br>
        !>  This error does not occur when the library is compiled with \ifx{2025.0.0 20241008}.<br>
        !>  It seems like this error occurs because of placing the following typed variable `speedup`
        !>  in the `forkjoin_parallelism_block` below.<br>
        !>  \remedy{2.0.0}
        !>  For now, the type definition and the typed variable declaration are taken out of the block and placed below.<br>
        !>  This must be checked with newer Intel compilers as `ifort` is being phased out by Intel.<br>
        type :: scaling_type
            real(RKG) :: maxval
            integer(IK) :: maxloc
            real(RKG), allocatable :: scaling(:)
            integer(IK), allocatable :: numproc(:)
        end type
        type :: speedup_type
            type(scaling_type) :: sameeff, zeroeff
        end type
        type(speedup_type) :: speedup

        character(*,SKG), parameter :: PROCEDURE_NAME = MODULE_NAME//SKG_"@getErrSampling()"
       !character(*,SKG), parameter :: NL1 = new_line(SKG_"a"), NL2 = NL1//NL1
        integer(IK) :: ndimp1, idim, iq
        type(proposal_type) :: proposal
        type(spec_type) :: spec
        type(stat_type) :: stat
        ndimp1 = ndim + 1_IK
        !!#if CAF_ENABLED
        !!! This must be coarray allocatable `[:]`, even though it will have only one element.
        !!! Otherwise, GNU Fortran 10 compiler results in runtime segmentation
        !!! fault in coarray mode, when the routine is called multiple times.
        !!! This used to be the case for a compound object. Does it also hold for intrinsic type objects?
        !!integer(IK), allocatable :: comv_unifRandState(:,:)[:]
        !!#else
        !!integer(IK), allocatable :: comv_unifRandState(:,:)
        !!#endif
        !
        !! Initialize the simulation auxiliary variables.
        !
        !!allocate(character(2**13 - 1, SK) :: errmsg) ! 8191: roughly 66Kb of memory for error message accumulation.
        !!reportFileUnit = output_unit ! Temporarily set the report file to stdout until the report file is set up.

        ! Setup the simulation specifications.

        spec = spec_type(modelr_type(0._RKG), method, ndim)
        err = spec%set(sampler)
        RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)

        !print *, """"//getStr([css_type(chainFileColName), css_type(spec%domainAxisName%val)], format = spec%chainFile%format%header)//achar(0, SKG)//""""
        if (spec%run%is%new .and. spec%image%is%leader) then
            if (spec%outputChainFileFormat%isBinary) then
                write(spec%chainFile%unit) getStr([css_type(chainFileColName), css_type(spec%domainAxisName%val)], format = spec%chainFile%format%header)//achar(0, SKG)
            else
                write(spec%chainFile%unit, spec%chainFile%format%header) (trim(chainFileColName(idim)), idim = 1, size(chainFileColName, 1, IK)), (trim(spec%domainAxisName%val(idim)), idim = 1, ndim)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        proposal = proposal_type(spec, err)
        RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (spec%image%is%leader) call spec%disp%text%wrap(NL1//spec%method%val//SKG_" simulation began on "//getDateTime(SK_"%c %z")//NL1)
        err = getErrKernelRun(getLogFunc, spec, stat, proposal)
        RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)
        if (spec%image%is%leader) call spec%disp%text%wrap(NL1//spec%method%val//SKG_" simulation ended on "//getDateTime(SK_"%c %z")//NL1)
        if (spec%image%is%first) then
            call spec%disp%skip(count = 2_IK, unit = output_unit)
            !call execute_command_line(" ")
            flush(output_unit)
        end if
        RETURN_IF_FAILED(__LINE__,FAILED_IMAGE(err%occurred),err%msg)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockLeaderPostProcessing: if (spec%image%is%leader) then

            ! Rules:
            !   -   Statistics values that require more than 1 line of the output file
            !       must be in the form of a table, the first line always being the column names.
            !       If the subsequent lines have one more column than the header line,
            !       then the first column of the subsequent lines will be interpreted as row names.
            !       Respecting this rule is important for parsing the contents of the report file in dynamic programming languages.

            associate(format => spec%reportFile%format%generic)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.numFuncCall.accepted")
                call spec%disp%show(stat%numFunCallAccepted, format = format)
                call spec%disp%note%show("This is the total number of accepted function calls (unique samples).")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.numFuncCall.acceptedRejected")
                call spec%disp%show(stat%numFunCallAcceptedRejected, format = format)
                call spec%disp%note%show("This is the total number of accepted or rejected function calls.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.numFuncCall.acceptedRejectedDelayed")
                call spec%disp%show(stat%numFunCallAcceptedRejectedDelayed, format = format)
                call spec%disp%note%show("This is the total number of accepted or rejected or delayed-rejection (if any requested) function calls.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.numFuncCall.acceptedRejectedDelayedUnused")
                call spec%disp%show(stat%numFunCallAcceptedRejectedDelayedUnused, format = format)
                call spec%disp%note%show("This is the total number of accepted or rejected or unused function calls (by all processes, including delayed rejections, if any requested).")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.efficiency.meanAcceptanceRate"))
                SET_DRAMDISE(call spec%disp%show(stat%cfc%meanAcceptanceRate(stat%numFunCallAccepted), format = format))
                SET_DRAMDISE(call spec%disp%note%show(SKG_"This is the average MCMC acceptance rate of the "//spec%method%val//SKG_" sampler."))

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.efficiency.acceptedOverAcceptedRejected"))
                SET_ParaNest(call spec%disp%show("stats.chain.uniques.efficiency.acceptedOverAcceptedRejected"))
                call spec%disp%show(real(stat%numFunCallAccepted, RKG) / real(stat%numFunCallAcceptedRejected, RKG), format = format) ! accepted2AcceptedRejected
                spec%msg = SKG_"This is the "//spec%method%val//SKG_" sampler efficiency given the accepted and rejected function calls, &
                &that is, the number of accepted function calls divided by the number of (accepted + rejected) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayed"))
                SET_ParaNest(call spec%disp%show("stats.chain.uniques.efficiency.acceptedOverAcceptedRejectedDelayed"))
                call spec%disp%show(real(stat%numFunCallAccepted, RKG) / real(stat%numFunCallAcceptedRejectedDelayed, RKG), format = format)
                spec%msg = SKG_"This is the "//spec%method%val//SKG_" sampler efficiency given the accepted, rejected, and delayed-rejection (if any requested) &
                &function calls, that is, the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayedUnused"))
                SET_ParaNest(call spec%disp%show("stats.chain.uniques.efficiency.acceptedOverAcceptedRejectedDelayedUnused"))
                call spec%disp%show(real(stat%numFunCallAccepted, RKG) / real(stat%numFunCallAcceptedRejectedDelayedUnused, RKG), format = format)
                spec%msg = SKG_"This is the "//spec%method%val//SKG_" sampler efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls &
                &(in parallel simulations), that is, the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.total")
                call spec%disp%show(stat%timer%clock, format = format)
                call spec%disp%note%show("This is the total runtime in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perFuncCallAccepted")
                call spec%disp%show(stat%timer%clock / real(stat%numFunCallAccepted, RKG), format = format)
                call spec%disp%note%show("This is the average effective time cost of each accepted function call, in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perFuncCallAcceptedRejected")
                call spec%disp%show(stat%timer%clock / real(stat%numFunCallAcceptedRejected, RKG), format = format)
                call spec%disp%note%show("This is the average effective time cost of each accepted or rejected function call, in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perFuncCallAcceptedRejectedDelayed")
                call spec%disp%show(stat%timer%clock / real(stat%numFunCallAcceptedRejectedDelayed, RKG), format = format)
                call spec%disp%note%show("This is the average effective time cost of each accepted or rejected function call (including delayed-rejections, if any requested), in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perFuncCallAcceptedRejectedDelayedUnused")
                call spec%disp%show(stat%timer%clock / real(stat%numFunCallAcceptedRejectedDelayedUnused, RKG), format = format)
                spec%msg = "This is the average effective time cost of each accepted or rejected or unused function call (including delayed-rejections, if any requested), in seconds. &
                &This timing can be directly compared to the corresponding timing of other parallel runs of the same simulation but with different processor counts to assess the parallel scaling. &
                &A lower timing value implies a better parallel scaling and performance."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perFuncCall")
                call spec%disp%show(stat%avgTimePerFunCall, format = format)
                call spec%disp%note%show("This is the average pure time cost of each function call, in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.time.perInterProcessCommunication")
                call spec%disp%show(stat%avgCommPerFunCall, format = format)
                call spec%disp%note%show("This is the average time cost of parallel inter-process communications per used (accepted or rejected or delayed-rejection) function call, in seconds.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin fork-join parallelism section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                forkjoin_parallelism_block: block

                    integer(IK) :: iell
                    integer(IK), parameter :: nscol = 5_IK

                    !!!!
                    !!!! First compute the fork-join strong scaling for the same and zero efficiency scenarios.
                    !!!!

                    ! \todo
                    ! The specified epsilon values for `parSecTime` and `comSecTime` are points of weakness of the following call
                    ! if the computers of the future civilization are roughly 1 billion times or more faster than the current 2020 technologies.
                    ! Therefore, a more robust solution is required for cases where the entire simulation is a dry run of the old existing simulation.
                    ! This situation is, however, such a rare occurrence that does not merit further investment in the current version of the library.
                    call setForkJoinScaling ( conProb = real(stat%numFunCallAccepted, RKG) / real(stat%numFunCallAcceptedRejected, RKG) & ! current sampling efficiency. ! LCOV_EXCL_LINE
                                            , seqSecTime = epsilon(1._RKG) & ! LCOV_EXCL_LINE time cost of the sequential section of the code, which is negligible here
                                            , comSecTime = real(stat%avgCommPerFunCall, RKG) / spec%image%count & ! LCOV_EXCL_LINE
                                            , parSecTime = real(stat%avgTimePerFunCall, RKG) & ! LCOV_EXCL_LINE
                                            , scalingMinLen = max(10_IK, spec%image%count * 2_IK) & ! LCOV_EXCL_LINE
                                            , scalingMaxLoc = speedup%sameeff%maxloc & ! LCOV_EXCL_LINE
                                            , scalingMaxVal = speedup%sameeff%maxval & ! LCOV_EXCL_LINE
                                            , numproc = speedup%sameeff%numproc & ! LCOV_EXCL_LINE
                                            , scaling = speedup%sameeff%scaling & ! LCOV_EXCL_LINE
                                            )
                    call setForkJoinScaling ( conProb = 0._RKG & ! zero sampling efficiency. ! LCOV_EXCL_LINE
                                            , seqSecTime = epsilon(1._RKG) & ! LCOV_EXCL_LINE time cost of the sequential section of the code, which is negligible here
                                            , comSecTime = real(stat%avgCommPerFunCall, RKG) / spec%image%count & ! LCOV_EXCL_LINE
                                            , parSecTime = real(stat%avgTimePerFunCall, RKG) & ! LCOV_EXCL_LINE
                                            , scalingMinLen = max(10_IK, spec%image%count * 2_IK) & ! LCOV_EXCL_LINE
                                            , scalingMaxLoc = speedup%zeroeff%maxloc & ! LCOV_EXCL_LINE
                                            , scalingMaxVal = speedup%zeroeff%maxval & ! LCOV_EXCL_LINE
                                            , numproc = speedup%zeroeff%numproc & ! LCOV_EXCL_LINE
                                            , scaling = speedup%zeroeff%scaling & ! LCOV_EXCL_LINE
                                            )

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.process.count.current")
                    call spec%disp%show(spec%image%count, format = format)
                    call spec%disp%note%show("This is the number of images/processes/threads used in this simulation.")

#if                 CAF_ENABLED || MPI_ENABLED || OMP_ENABLED
                    if (spec%image%count == 1_IK) then
                        spec%msg = spec%method%val//SKG_" is used in parallel mode with only one image/process/thread. This can be computationally inefficient. &
                        &Consider using the serial version of the code or provide more processes at runtime if it seems to be beneficial as discussed below."
                        call spec%disp%note%show(spec%msg, unit = output_unit)!, format = format, tmsize = 3_IK, bmsize = 1_IK)
                    end if
#endif

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.process.count.optimal.sameeff")
                    call spec%disp%show(speedup%sameeff%maxloc, format = format)
                    spec%msg = SKG_"This is the predicted optimal number of physical computing processes for "//spec%parallelism%val// & ! LCOV_EXCL_LINE
                    SKG_" parallelization model, assuming the same (current) sampling efficiency and parallel communication overhead as in this simulation."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.process.count.optimal.zeroeff")
                    call spec%disp%show(speedup%zeroeff%maxloc, format = format)
                    spec%msg = "This is the predicted number of physical computing processes for "//spec%parallelism%val// & ! LCOV_EXCL_LINE
                    SKG_" parallelization model, assuming zero sampling efficiency and the same (current) parallel communication overhead as in this simulation. &
                    &This sampling task will likely NOT benefit from any additional computing processes beyond this predicted optimal count, "// & ! LCOV_EXCL_LINE
                    getStr(speedup%zeroeff%maxloc)//SKG_", in the above, under the ideal synchronous fork-join parallelism scenario. &
                    &This is true for any value of sampling efficiency given the current parallel communication overhead. &
                    &Keep in mind that the predicted optimal number of processes in this zero-efficiency sampling scenario is only an &
                    &estimate whose accuracy depends on many runtime factors, including the topology of the parallel communication network used, &
                    &the number of processes per node, and the number of tasks to each processor or node."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.speedup.current")
                    call spec%disp%show(speedup%sameeff%scaling(spec%image%count), format = format)
                    spec%msg = "This is the estimated maximum speedup gained via "//spec%parallelism%val// & ! LCOV_EXCL_LINE
                    SKG_" parallelization model compared to serial mode given the current parallel communication overhead."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    spec%msg = ""
                    call spec%disp%show("stats.parallelism.speedup.optimal.sameeff")
                    call spec%disp%show(speedup%sameeff%maxval, format = format)
                    if (spec%parallelism%is%forkJoin .and. speedup%sameeff%scaling(spec%image%count) < 1._RKG) then
                        spec%msg = "The time cost of calling the user-provided objective function must be at least "//getStr(1._RKG / speedup%sameeff%scaling(spec%image%count), SK_"(g0.6)")//&
                        SKG_" times more (that is, ~"//getStr(10**6 * stat%avgTimePerFunCall / speedup%sameeff%scaling(spec%image%count), SK_"(g0.6)")//" microseconds) to see any performance benefits from "//&
                        spec%parallelism%val//SKG_" parallelization model for this simulation. "
                        if (speedup%sameeff%maxloc == 1_IK) then
                            spec%msg = spec%msg//SKG_"Consider simulating in serial mode in the future (if used within the same computing environment and with the same configuration as used here)."
                        else
                            spec%msg = spec%msg//SKG_"Consider simulating on "//getStr(speedup%sameeff%maxloc)//&
                            SKG_" processes in the future (if used within the same computing environment and with the same configuration as used here)."
                        end if
                        if (.not. spec%outputSplashMode%is%silent) call spec%disp%note%show(spec%msg, unit = output_unit, tmsize = spec%disp%tmsize, bmsize = spec%disp%bmsize)
                        spec%msg = NL1//spec%msg
                    end if
                    spec%msg = SKG_"This is the predicted optimal maximum speedup gained via "//spec%parallelism%val// & ! LCOV_EXCL_LINE
                    SKG_" parallelization model, given the current sampling efficiency and parallel communication overhead."//spec%msg
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.speedup.optimal.zeroeff")
                    call spec%disp%show(speedup%zeroeff%maxval, format = format)
                    spec%msg = SKG_"This is the predicted optimal maximum speedup gained via `"//spec%parallelism%val// & ! LCOV_EXCL_LINE
                    SKG_"` parallelization model, assuming the ideal zero-efficiency sampling and &
                    &the same (current) parallel communication overhead as in this simulation."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    !!!!
                    !!!! Find individual image contributions and the Cyclic Geometric fit to the contributions.
                    !!!!

                    imageContribution_block: block

                        integer(IK) :: pidSuccessLen
                        real(RKG) :: successProbNormFac(2)
                        integer(IK), allocatable :: index(:), cntSuccess(:), pidSuccess(:)

                        call setResized(index, spec%image%count)
                        call setResized(cntSuccess, spec%image%count)
                        call setResized(pidSuccess, spec%image%count)
                        call setUnique(stat%cfc%processID, unique = pidSuccess, lenUnique = pidSuccessLen, count = cntSuccess)

                        if (pidSuccessLen < spec%image%count) then
                            pidSuccess(pidSuccessLen + 1 :) = getComplementRange(pidSuccess(1 : pidSuccessLen), start = 1_IK, stop = spec%image%count, step = 1_IK)
                            cntSuccess(pidSuccessLen + 1 :) = 0_IK
                        end if
                        call setSorted(pidSuccess, index)
                        pidSuccess(:) = pidSuccess(index)
                        cntSuccess(:) = cntSuccess(index)

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        call spec%disp%show("stats.parallelism.process.contribution.count.current")
                        call spec%disp%show(css_type([character(15) :: "processID", "count"]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                        do iell = 1, spec%image%count
                            write(spec%reportFile%unit, spec%reportFile%format%integer) iell, cntSuccess(iell)
                        end do
                        call spec%disp%skip(count = spec%disp%bmsize)
                        spec%msg = SKG_"These are the contributions of individual processes to the construction of the output chain of the "//spec%method%val//SKG_" sampler. &
                        &Each count value represents the total number of accepted states (useful contributions to the simulation) by the corresponding processor, starting &
                        &from the first processor to the last. This information is mainly useful in synchronous parallel Fork-Join (singleChain) simulations. &
                        &Ideally, one would desire equal contributions from all processes to the final output chain, although this is never the case."
                        call spec%disp%note%show(spec%msg)

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        call spec%disp%show("stats.parallelism.process.contribution.count.current.fit")
                        call spec%disp%show(css_type([character(15) :: "successProb", "normFac"]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                        if (spec%parallelism%is%forkJoin) then
                            pidSuccess = pack(pidSuccess, 0 < cntSuccess)
                            cntSuccess = pack(cntSuccess, 0 < cntSuccess)
                            err%occurred = isFailedGeomCyclicFit(pidSuccess, cntSuccess, spec%image%count, successProbNormFac(1), successProbNormFac(2))
                            if (err%occurred) then
                                call spec%disp%show(css_type([character(15) :: "NaN", "NaN"]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                                err%occurred = .false._LK
                            else
                                call spec%disp%show(successProbNormFac, format = spec%reportFile%format%allreal, tmsize = 0_IK)
                            end if
                        else
                            successProbNormFac = [1._RKG, real(cntSuccess(size(cntSuccess, 1, IK)), RKG)]
                            call spec%disp%show(successProbNormFac, format = spec%reportFile%format%allreal, tmsize = 0_IK)
                        end if
                        call spec%disp%skip(count = spec%disp%bmsize)
                        spec%msg =  SKG_"These are the two parameters of the Cyclic Geometric fit to the distribution of the processor contributions to the construction &
                                        &of the final output chain of visited states. The processor contributions are reported in the first column of the output chain file. &
                                        &The processor contribution frequencies are listed above. The fit has the following form: "//NL2// &
                                    SKG_"    processConstribution(i) ="//NL1//&
                                        SKG_"successProb * normFac * (1 - successProb)^(i - 1) / (1 - (1 - successProb)^processCount)"//NL2// &
                                    SKG_"where `i` is the ID of the processor (starting from index `1`), `processCount` is the total number of &
                                        &processes used in the simulation, `successProb` is equivalent to an effective sampling efficiency computed &
                                        &from the contributions of individual processes to the output chain, and `normFac` is a normalization constant."
                        call spec%disp%note%show(spec%msg)

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    end block imageContribution_block

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.speedup.scaling.strong.sameeff")
                    call spec%disp%show(css_type([character(15) :: "processCount", "speedup"]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    do iell = 1, size(speedup%sameeff%scaling, 1, IK)
                        write(spec%reportFile%unit, spec%reportFile%format%intreal) speedup%sameeff%numproc(iell), speedup%sameeff%scaling(iell)
                        !call spec%disp%show(speedup%sameeff%scaling(iell : min(iell + nscol - 1_IK, size(speedup%sameeff%scaling, 1, IK))), format = scalingFormat, tmsize = 0_IK, bmsize = 0_IK)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                    spec%msg = SKG_"This is the predicted strong-scaling speedup behavior of the "//spec%parallelism%val//SKG_" parallelization model, &
                    &given the current parallel communication overhead in the above and the current sampling efficiency, for increasing numbers of processes, starting from a single process."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.parallelism.speedup.scaling.strong.zeroeff")
                    call spec%disp%show(css_type([character(15) :: "processCount", "speedup"]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    do iell = 1, size(speedup%zeroeff%scaling, 1, IK)
                        write(spec%reportFile%unit, spec%reportFile%format%intreal) speedup%zeroeff%numproc(iell), speedup%zeroeff%scaling(iell)
                        !call spec%disp%show(speedup%zeroeff%scaling(iell : min(iell + nscol - 1_IK, size(speedup%zeroeff%scaling, 1, IK))), format = scalingFormat, tmsize = 0_IK, bmsize = 0_IK)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                    spec%msg = SKG_"This is the predicted strong-scaling speedup behavior of the "//spec%parallelism%val//SKG_" parallelization model, &
                    &assuming an ideal synchronous fork-join parallelism scenario (with 0% sampling efficiency), given the current parallel communication overhead, &
                    &for increasing numbers of processes, starting from a single process."
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end block forkjoin_parallelism_block

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.domain.logVolume")
                call spec%disp%show(spec%domain%logVol, format = format)
                call spec%disp%note%show("This is the natural logarithm of the volume of the "//spec%ndim%str//SKG_"-dimensional domain over which the objective function was defined.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if             ParaDRAM_ENABLED || ParaDISE_ENABLED
                if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                    call spec%disp%note%show("Computing the statistical properties of the Markov chain...", unit = output_unit)
                end if
                call spec%disp%text%wrap(NL1//SKG_"The statistical properties of the Markov chain"//NL1)
#elif           ParaNest_ENABLED
                if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                    call spec%disp%note%show("Computing the statistical properties of the output chain...", unit = output_unit)
                end if
                call spec%disp%text%wrap(NL1//SKG_"The statistical properties of the output chain"//NL1)
#else
#error          "Unrecognized sampler."
#endif
                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.logFunc.max.val"))
                SET_ParaNest(call spec%disp%show("stats.chain.uniques.logFunc.max.val"))
                call spec%disp%show(stat%chain%mode%val, format = format)
                call spec%disp%note%show("This is the maximum logFunc value (the maximum of the user-specified objective function).")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SET_DRAMDISE(call spec%disp%show("stats.chain.verbose.logFunc.max.crd"))
                SET_ParaNest(call spec%disp%show("stats.chain.uniques.logFunc.max.crd"))
                call spec%disp%show(css_type(spec%domainAxisName%val), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%show(stat%chain%mode%crd, format = spec%reportFile%format%allreal)
                call spec%disp%note%show("This is the coordinates, within the domain of the user-specified objective function, where the maximum function value occurs.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if             ParaDRAM_ENABLED || ParaDISE_ENABLED
                call spec%disp%show("stats.chain.compact.logFunc.max.loc")
                call spec%disp%show(stat%chain%mode%loc, format = format)
                call spec%disp%note%show("This is the location of the first occurrence of the maximum logFunc in the compact chain.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.logFunc.max.loc")
                call spec%disp%show(sum(stat%cfc%sampleWeight(1 : stat%chain%mode%loc - 1)) + 1_IK, format = format) ! stat%chain%mode%Loc%verbose
                call spec%disp%note%show("This is the location of the first occurrence of the maximum logFunc in the verbose (Markov) chain.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.compact.burnin.loc.likelihoodBased")
                call spec%disp%show(stat%burninLocMCMC%compact, format = format)
                call spec%disp%note%show("This is the burnin location in the compact chain, based on the occurrence likelihood.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.compact.burnin.loc.adaptationBased")
                call spec%disp%show(stat%burninLocDRAM%compact, format = format)
                call spec%disp%note%show("This is the burnin location in the compact chain, based on the value of proposalAdaptationBurnin simulation specification.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.burnin.loc.likelihoodBased")
                call spec%disp%show(stat%burninLocMCMC%verbose, format = format)
                call spec%disp%note%show("This is the burnin location in the verbose (Markov) chain, based on the occurrence likelihood.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.burnin.loc.adaptationBased")
                call spec%disp%show(stat%burninLocDRAM%verbose, format = format)
                call spec%disp%note%show("This is the burnin location in the verbose (Markov) chain, based on the value of proposalAdaptationBurnin simulation specification.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! reset burninLocation to the maximum value

                if (stat%burninLocDRAM%compact > stat%burninLocMCMC%compact) then
                    stat%burninLocMCMC%compact = stat%burninLocDRAM%compact
                    stat%burninLocMCMC%verbose = stat%burninLocDRAM%verbose
                end if

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Compute the statistical properties of the MCMC chain
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! Compute the covariance and correlation upper-triangle matrices.

                !>  \warning
                !>  forrtl: severe (174): SIGSEGV, segmentation fault occurred
                !>  Apparently, when the chain is too long (e.g., 300000), the stack size can lead to such an error.
                !>  The stack size will have to be adjusted in such cases.
                !>  \todo
                !>  A robust fix to the above segmentation fault error required further digging into this problem.
                !>  This could become a serious hard-to-resolve error at higher levels in dynamic languages, especially on Windows.
                !>  update: The problem is now resolved by allocating on the heap.

                call setResized(stat%chain%avg, ndim)
                call setResized(stat%chain%cov, [ndim, ndim])
                call setMean(stat%chain%avg, stat%cfc%sampleState(:, stat%burninLocMCMC%compact : stat%cfc%nsam), 2_IK, stat%cfc%sampleWeight(stat%burninLocMCMC%compact : stat%cfc%nsam), stat%chain%size)
                call setCov(stat%chain%cov, lowDia, stat%chain%avg, stat%cfc%sampleState(:, stat%burninLocMCMC%compact : stat%cfc%nsam), 2_IK, stat%cfc%sampleWeight(stat%burninLocMCMC%compact : stat%cfc%nsam), stat%chain%size)
                call setMatCopy(stat%chain%cov, rdpack, stat%chain%cov, rdpack, lowDia, transHerm)
                stat%chain%cor = getCor(stat%chain%cov, subsetv = lowDia)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.size.burninExcluded")
                call spec%disp%show(stat%chain%size, format = format)
                call spec%disp%note%show("This is the length of the verbose (Markov) Chain excluding burnin.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.avg")
                call spec%disp%show(css_type(spec%domainAxisName%val), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%show(stat%chain%avg, format = spec%reportFile%format%allreal)
                call spec%disp%note%show("This is the mean (avg) of the verbose (Markov) chain variables.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.std")
                call spec%disp%show(css_type(spec%domainAxisName%val), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%show([(sqrt(stat%chain%cov(idim, idim)), idim = 1, ndim)], format = spec%reportFile%format%allreal)
                call spec%disp%note%show("This is the standard deviation (std) of the verbose (Markov) chain variables.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.covmat")
                call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: SKG_"axis", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%skip(count = spec%disp%tmsize)
                do idim = 1, ndim
                    write(spec%reportFile%unit, spec%reportFile%format%strreal) trim(adjustl(spec%domainAxisName%val(idim))), stat%chain%cov(1 : ndim, idim)
                end do
                call spec%disp%skip(count = spec%disp%bmsize)
                call spec%disp%note%show("This is the covariance matrix of the verbose (Markov) chain.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.verbose.cormat")
                call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: SKG_"axis", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%skip(count = spec%disp%tmsize)
                do idim = 1, ndim
                    write(spec%reportFile%unit, spec%reportFile%format%strreal) trim(adjustl(spec%domainAxisName%val(idim))), stat%chain%cor(1 : ndim, idim)
                end do
                call spec%disp%skip(count = spec%disp%bmsize)
                call spec%disp%note%show("This is the correlation matrix of the verbose (Markov) chain.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! Compute the chain quantiles.
                ! \bug Intel ifort bug with heap memory allocations:
                ! The Intel ifort cannot digest the following task for an input 2D sample, throwing `double free or corruption (out)`.
                ! The source of this error was traced back to the return point from `setExtrap()` within `getQuan()`.
                ! The runtime error seems to resolve if the quantile matrix is computed in a loop like the following.

                call setResized(stat%chain%quantile%quan, [size(stat%chain%quantile%prob, 1, IK), ndim])
                block
                    real(RKG), allocatable :: sample(:)
                    do idim = 1, ndim
                        sample = stat%cfc%sampleState(idim, stat%burninLocMCMC%compact : stat%cfc%nsam)
                        stat%chain%quantile%quan(:, idim) = getQuan(neimean, stat%chain%quantile%prob, sample, stat%cfc%sampleWeight(stat%burninLocMCMC%compact : stat%cfc%nsam), stat%chain%size)
                    end do
                end block
                !stat%chain%quantile%quan = getQuan(neimean, stat%chain%quantile%prob, stat%cfc%sampleState(idim, stat%burninLocMCMC%compact : stat%cfc%nsam), 2_IK, stat%cfc%sampleWeight(stat%burninLocMCMC%compact : stat%cfc%nsam), stat%chain%size)

                call spec%disp%show("stats.chain.verbose.quantile")
                call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: SKG_"quantile", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                call spec%disp%skip(count = spec%disp%tmsize)
                do iq = 1, size(stat%chain%quantile%prob, 1, IK)
                    write(spec%reportFile%unit, spec%reportFile%format%allreal) stat%chain%quantile%prob(iq), (stat%chain%quantile%quan(iq, idim), idim = 1, ndim)
                end do
                call spec%disp%skip(count = spec%disp%bmsize)
                call spec%disp%note%show("This is the quantiles table of the variables of the verbose (Markov) chain.")
#endif
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Generate the i.i.d. sample statistics and output file (if requested)
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! report refined sample statistics, and generate output refined sample if requested.

                if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                    call spec%disp%note%show("Computing the final refined i.i.d. sample size...", unit = output_unit)
                end if
#if             ParaDRAM_ENABLED || ParaDISE_ENABLED
                err = stat%sfc%getErrRefinement ( sampleWeight = stat%cfc%sampleWeight(stat%burninLocMCMC%compact :) & ! LCOV_EXCL_LINE
                                                , sampleLogFunc = stat%cfc%sampleLogFunc(stat%burninLocMCMC%compact :) & ! LCOV_EXCL_LINE
                                                , sampleState = stat%cfc%sampleState(:, stat%burninLocMCMC%compact :) & ! LCOV_EXCL_LINE
                                                , outputSampleRefinementCount = spec%outputSampleRefinementCount%val & ! LCOV_EXCL_LINE
                                                , outputSampleRefinementMethod = spec%outputSampleRefinementMethod%val & ! LCOV_EXCL_LINE
                                                )

                ! compute the maximum integrated autocorrelation times for each variable.

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.act")
                call spec%disp%skip(count = spec%disp%tmsize)
                write(spec%reportFile%unit, spec%reportFile%format%fixform) "refinementStage", "sampleSize", "ACT_sampleLogFunc", (SKG_"ACT_"//trim(adjustl(spec%domainAxisName%val(idim))), idim = 1, ndim)
                block
                    integer(IK) :: iref
                    character(:,SKG), allocatable :: thisForm
                    thisForm = SKG_"('"//spec%reportFile%indent//SKG_"',2(I"//spec%outputColumnWidth%max//SKG_",' '),*("//& ! LCOV_EXCL_LINE
                    spec%real%ed//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//spec%real%ex//SKG_",:,' '))"
                    call spec%disp%skip(count = spec%disp%tmsize)
                    do iref = 0, stat%sfc%nref
                        write(spec%reportFile%unit, thisForm) iref, stat%sfc%sumw(iref), stat%sfc%act(0 : ndim, iref)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                end block
                spec%msg = SKG_"This is the table of the Integrated Autocorrelation (ACT) of individual variables in the verbose (Markov) chain, at increasing stages of chain refinements."
                if (stat%sfc%nref == 0_IK) spec%msg = spec%msg//NL1// & ! LCOV_EXCL_LINE
                SKG_"The user-specified `outputSampleRefinementCount` ("//getStr(spec%outputSampleRefinementCount%val)//SKG_") &
                &is too small to ensure an accurate computation of the decorrelated i.i.d. effective sample size. No refinement of the Markov chain was performed."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! Report the final Effective Sample Size (ESS) based on ACT.

                stat%ess = stat%sfc%sumw(stat%sfc%nref)

#elif           ParaNest_ENABLED

                ! Read the ParaNest output chain file contents.

                err = getErrChainRead(stat%cfc, file = spec%chainFile%file, spec%outputChainFileFormat%val, pre = stat%numFunCallAccepted)
                if (err%occurred) then
                    err%msg = PROCEDURE_NAME//err%msg ! LCOV_EXCL_LINE
                    exit blockLeaderPostProcessing ! LCOV_EXCL_LINE
                    return ! LCOV_EXCL_LINE
                end if

                ! Compute the effective sample size.

                stat%ess = nint(sum(exp(stat%cfc%sampleLogWeight - maxval(stat%cfc%sampleLogWeight))), IK)
#else
#error          "Unrecognized interface."
#endif
                stat%sample%size = abs(spec%outputSampleSize%val)
                if (spec%outputSampleSize%val < 0_IK) stat%sample%size = stat%sample%size * stat%ess

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.ess")
                call spec%disp%show(stat%ess, format = format)
                call spec%disp%note%show("This is the estimated Effective (i.i.d.) Sample Size (ESS) of the final refined chain.")

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.efficiency.essOverAccepted")
                call spec%disp%show(real(stat%ess, RKG) / real(stat%numFunCallAccepted, RKG), format = format)
                spec%msg = SKG_"This is the effective sampling efficiency given the accepted function calls, that is, &
                &the final refined effective sample size (ESS) divided by the number of accepted function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.efficiency.essOverAcceptedRejected")
                call spec%disp%show(real(stat%ess, RKG) / real(stat%numFunCallAcceptedRejected, RKG), format = format)
                spec%msg = SKG_"This is the effective sampling efficiency given the accepted and rejected function calls, that is, &
                &the final refined effective sample size (ESS) divided by the number of (accepted + rejected) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.efficiency.essOverAcceptedRejectedDelayed")
                call spec%disp%show(real(stat%ess, RKG) / real(stat%numFunCallAcceptedRejectedDelayed, RKG), format = format)
                spec%msg = SKG_"This is the effective sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, &
                &that is, the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                call spec%disp%show("stats.chain.refined.efficiency.essOverAcceptedRejectedDelayedUnused")
                call spec%disp%show(real(stat%ess, RKG) / real(stat%numFunCallAcceptedRejectedDelayedUnused, RKG), format = format)
                spec%msg = SKG_"This is the effective sampling efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls, &
                &(in parallel simulations), that is, the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
                call spec%disp%note%show(spec%msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                !end associate blockEffectiveSampleSize

                ! Generate the output refined sample if requested.

                blockSampleFileGeneration: if (spec%outputSampleSize%val == 0_IK) then

                    call spec%disp%note%show(SKG_"Skipping the generation of the refined chain and the output "//spec%sampleFile%kind//SKG_" file, as requested by the user...")

                else blockSampleFileGeneration

                    ! report to the report-file(s)

                    call spec%disp%note%show(SKG_"Generating the output "//spec%sampleFile%kind//SKG_" file:"//NL1//spec%sampleFile%file)
                    if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                        ! print the message for the generating the output sample file on the first image
                        call spec%disp%note%show(SKG_"Generating the output "//spec%sampleFile%kind//SKG_" file:"//NL1//spec%sampleFile%file, unit = output_unit, bmsize = 0_IK)
                        ! print the message for the generating the output sample file on the rest of the images in order.
#if                     CAF_ENABLED || MPI_ENABLED || OMP_ENABLED
                        if (spec%parallelism%is%multiChain) then
                            block
                                integer(IK) :: imageID
                                do imageID = 2, spec%image%count
                                    call spec%disp%note%show(getReplaced(spec%sampleFile%file, SKG_"pid1", SKG_"pid"//getStr(imageID)), unit = output_unit, tmsize = 0_IK, bmsize = 0_IK)
                                end do
                            end block
                        end if
#endif
                        call spec%disp%skip(unit = output_unit)
                    end if

                    ! Begin sample file generation.

                    stat%sfc%colname = [css_type(chainFileColName(size(chainFileColName, 1, IK))), css_type(spec%domainAxisName%val)]
                    stat%sfc%header = getStr(stat%sfc%colname, format = spec%sampleFile%format%header)
#if                 ParaDISE_ENABLED || ParaDRAM_ENABLED
                    if (spec%outputSampleSize%val /= -1_IK) then
                        ! Regenerate the refined sample, this time with the user-specified sample size.
                        block
                            integer(IK), allocatable :: cumSumWeight(:), unifrnd(:)
                            integer(IK) :: isam, iloc
                            stat%sfc%nref = 0_IK
                            call setRebound(stat%sfc%nsam, 0_IK, 0_IK)
                            call setRebound(stat%sfc%sumw, 0_IK, 0_IK)
                            stat%sfc%nsam(stat%sfc%nref) = stat%numFunCallAccepted - stat%burninLocMCMC%compact + 1_IK
                            call setResized(cumSumWeight, stat%sfc%nsam(stat%sfc%nref))
                            call setCumSum(cumSumWeight, stat%cfc%sampleWeight(stat%burninLocMCMC%compact :))
                            unifrnd = getUnifRand(spec%rng, 1_IK, cumSumWeight(stat%sfc%nsam(stat%sfc%nref)), stat%sample%size)
                            call setRebound(stat%sfc%sampleLogFuncState, [0_IK, 1_IK], [ndim, stat%sfc%nsam(stat%sfc%nref)])
                            stat%sfc%sampleLogFuncState(1 :, :) = stat%cfc%sampleState(:, stat%burninLocMCMC%compact :)
                            stat%sfc%sampleLogFuncState(0, :) = stat%cfc%sampleLogFunc(stat%burninLocMCMC%compact :)
                            stat%sfc%sampleWeight = getFilled(0_IK, stat%sfc%nsam(stat%sfc%nref))
                            loopOverSample: do isam = 1, stat%sample%size
                                do iloc = 1, size(cumSumWeight, 1, IK)
                                    if (cumSumWeight(iloc) < unifrnd(isam)) cycle
                                    stat%sfc%sampleWeight(iloc) = stat%sfc%sampleWeight(iloc) + 1_IK
                                    cycle loopOverSample
                                end do
                                error stop "This is a miracle! Please report it to the developers at: https://github.com/cdslaborg/paramonte"
                            end do loopOverSample
                            stat%sfc%sumw(stat%sfc%nref) = sum(stat%sfc%sampleWeight)
                        end block
                    end if
#elif               ParaNest_ENABLED
                    stat%sfc%nref = 0_IK
                    call setRebound(stat%sfc%nsam, 0_IK, 0_IK)
                    call setRebound(stat%sfc%sumw, 0_IK, 0_IK)
                    stat%sfc%sampleWeight = nint(getReweight(exp(stat%cfc%sampleLogWeight - getLogSumExp(stat%cfc%sampleLogWeight)), 1._RKG / stat%sample%size), IK)
                    stat%sfc%nsam(stat%sfc%nref) = size(stat%sfc%sampleWeight, 1, IK)
                    call setRebound(stat%sfc%sampleLogFuncState, [0_IK, 1_IK], [ndim, stat%sfc%nsam(stat%sfc%nref)])
                    stat%sfc%sumw(stat%sfc%nref) = sum(stat%sfc%sampleWeight)
#else
#error              "Unrecognized interface."
#endif
                    ! open the output sample file and write the sample.

                    sfc_block: block
                        integer(IK) :: isam, iwei
                        spec%sampleFile%iomsg = SKG_""
                        spec%sampleFile%unit = getFileUnit() ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
                        open(unit = spec%sampleFile%unit, file = spec%sampleFile%file, form = spec%sampleFile%form, status = spec%sampleFile%status, position = spec%sampleFile%position, iostat = spec%sampleFile%iostat, iomsg = spec%sampleFile%iomsg SHARED)
                        if (spec%sampleFile%iostat /= 0_IK) exit sfc_block
                        write(spec%sampleFile%unit, spec%sampleFile%format%header, iostat = spec%sampleFile%iostat, iomsg = spec%sampleFile%iomsg) stat%sfc%header
                        if (spec%sampleFile%iostat /= 0_IK) exit sfc_block
                        do isam = 1, stat%sfc%nsam(stat%sfc%nref)
                            do iwei = 1, stat%sfc%sampleWeight(isam)
                                write(spec%sampleFile%unit, spec%sampleFile%format%rows, iostat = spec%sampleFile%iostat, iomsg = spec%sampleFile%iomsg) stat%sfc%sampleLogFuncState(:, isam)
                                if (spec%sampleFile%iostat /= 0_IK) exit sfc_block
                            end do
                        end do
                        close(spec%sampleFile%unit, iostat = spec%sampleFile%iostat, iomsg = spec%sampleFile%iomsg)
                        if (spec%sampleFile%iostat /= 0_IK) exit sfc_block
                    end block sfc_block
                    err%occurred = spec%sampleFile%iostat /= 0_IK
                    if (err%occurred) then
                        err%stat = spec%sampleFile%iostat
                        err%msg = spec%sampleFile%iomsg
                        exit blockLeaderPostProcessing
                    end if

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Compute the statistical properties of the refined sample
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                        call spec%disp%note%show(SKG_"Computing the statistical properties of the final output sample...", unit = output_unit)
                    end if
                    call spec%disp%text%wrap(NL1//SKG_"The statistical properties of the final output sample"//NL1)
#if                 ParaDRAM_ENABLED || ParaDISE_ENABLED
                    CHECK_ASSERTION(__LINE__, stat%sample%size == stat%sfc%sumw(stat%sfc%nref), \
                    SK_"@getErrSampling(): The condition `stat%sample%size == stat%sfc%sumw(stat%sfc%nref)` must hold. stat%sample%size, stat%sfc%sumw(stat%sfc%nref), stat%sfc%nref = "//\
                    getStr([stat%sample%size, stat%sfc%sumw(stat%sfc%nref), stat%sfc%nref]))
#endif
                    ! Compute the covariance and correlation upper-triangle matrices.

                    call setRebound(stat%sample%avg, 0_IK, ndim)
                    call setRebound(stat%sample%cov, [0_IK, 0_IK], [ndim, ndim])
                    call setMean(stat%sample%avg, stat%sfc%sampleLogFuncState(:, 1 : stat%sfc%nsam(stat%sfc%nref)), 2_IK, stat%sfc%sampleWeight(1 : stat%sfc%nsam(stat%sfc%nref)), stat%sfc%sumw(stat%sfc%nref))
                    if (ndim < stat%sfc%nsam(stat%sfc%nref)) then
                        call setCov(stat%sample%cov, lowDia, stat%sample%avg, stat%sfc%sampleLogFuncState(:, 1 : stat%sfc%nsam(stat%sfc%nref)), 2_IK, stat%sfc%sampleWeight(1:stat%sfc%nsam(stat%sfc%nref)), stat%sfc%sumw(stat%sfc%nref))
                        call setMatCopy(stat%sample%cov, rdpack, stat%sample%cov, rdpack, lowDia, transHerm)
                        stat%sample%cor = getCor(stat%sample%cov, subsetv = lowDia)
                    else
                        stat%sample%cor = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 0._RKG)
                        stat%sample%cov = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 0._RKG)
                    end if

                    ! Report the refined chain statistics.

                    !formatStr = "(1A"//spec%outputColumnWidth%max//SKG_",*(E"//spec%outputColumnWidth%max//SKG_"."//spec%outputPrecision%str//SKG_"))"

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.chain.refined.length")
                    call spec%disp%show(stat%sample%size, format = format)
                    spec%msg = "This is the final output refined sample size."
                    if (spec%outputSampleSize%val /= -1_IK) then
                        if (stat%sample%size < stat%ess) then
                            spec%msg = spec%msg//NL1//SKG_"The user-requested sample size ("//getStr(stat%sample%size)//SKG_") is smaller &
                            &than the potentially-optimal i.i.d. sample size ("//getStr(stat%ess)//SKG_"). The output sample &
                            &contains i.i.d. samples, however, the sample-size could have been larger if it had been set to the optimal size. &
                            &To get the optimal size in the future runs, set outputSampleSize = -1, or drop it from the input list."
                        elseif (stat%ess < stat%sample%size) then
                            spec%msg = spec%msg//NL1//SKG_"The user-requested sample size ("//getStr(stat%sample%size)//SKG_") is larger than &
                            &the potentially-optimal i.i.d. sample size ("//getStr(stat%ess)//SKG_"). The resulting sample likely contains &
                            &duplicates and is not independently and identically distributed (i.i.d.). To get the optimal &
                            &size in the future runs, set outputSampleSize = -1, or drop it from the input list."
                        else ! LCOV_EXCL_LINE
                            spec%msg = spec%msg//NL1//SKG_"How lucky that could be! The user-requested sample size ("//getStr(stat%sample%size)//& ! LCOV_EXCL_LINE
                            SKG_") is equal to the potentially-optimal i.i.d. sample size determined by the "//spec%method%val//SKG_" sampler." ! LCOV_EXCL_LINE
                        end if
                    end if
                    call spec%disp%note%show(spec%msg)

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.chain.refined.avg")
                    call spec%disp%show(css_type(spec%domainAxisName%val), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    call spec%disp%show(stat%sample%avg(1:), format = spec%reportFile%format%allreal)
                    call spec%disp%note%show("This is the mean (avg) of the final output refined sample.")

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.chain.refined.std")
                    call spec%disp%show(css_type(spec%domainAxisName%val), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    call spec%disp%show([(sqrt(stat%sample%cov(idim, idim)), idim = 1, ndim)], format = spec%reportFile%format%allreal)
                    call spec%disp%note%show("This is the standard deviation (std) of the final output refined sample.")

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.chain.refined.covmat")
                    call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: SKG_"axis", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    call spec%disp%skip(count = spec%disp%tmsize)
                    do idim = 1, ndim
                        write(spec%reportFile%unit, spec%reportFile%format%strreal) trim(adjustl(spec%domainAxisName%val(idim))), stat%sample%cov(1 : ndim, idim)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                    call spec%disp%note%show("This is the covariance matrix of the final output refined sample.")
                    if (.not. isMatClass(stat%sample%cov, posdefmat)) call spec%disp%warn%show("The sample covariance matrix is not positive-definite. sample size = "//getStr(stat%sfc%nsam(stat%sfc%nref)))

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    call spec%disp%show("stats.chain.refined.cormat")
                    call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: SKG_"axis", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    call spec%disp%skip(count = spec%disp%tmsize)
                    do idim = 1, ndim
                        write(spec%reportFile%unit, spec%reportFile%format%strreal) trim(adjustl(spec%domainAxisName%val(idim))), stat%sample%cor(1 : ndim, idim)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                    call spec%disp%note%show("This is the correlation matrix of the final output refined sample.")
                    if (.not. isMatClass(stat%sample%cor, posdefmat)) call spec%disp%warn%show("The sample correlation matrix is not positive-definite. sample size = "//getStr(stat%sfc%nsam(stat%sfc%nref)))

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! Compute the sample quantiles.
                    ! \bug Intel ifort bug with heap memory allocations:
                    ! The Intel ifort cannot digest the following task for an input 2D sample, throwing `double free or corruption (out)`.
                    ! The source of this error was traced back to the return point from `setExtrap()` within `getQuan()`.
                    ! The runtime error seems to resolve if the quantile matrix is computed in a loop like the following.

                    call setResized(stat%sample%quantile%quan, [size(stat%sample%quantile%prob, 1, IK), ndim])
                    block
                        real(RKG), allocatable :: sample(:)
                        do idim = 1, ndim
                            sample = stat%sfc%sampleLogFuncState(idim, :)
                            stat%sample%quantile%quan(:, idim) = getQuan(neimean, stat%sample%quantile%prob, sample, stat%sfc%sampleWeight, stat%sfc%sumw(stat%sfc%nref))
                        end do
                    end block
                    !stat%sample%quantile%quan = getQuan(neimean, stat%sample%quantile%prob, stat%sfc%sampleLogFuncState, 2_IK, stat%sfc%sampleWeight, stat%sfc%sumw(stat%sfc%nref))

                    call spec%disp%show("stats.chain.refined.quantile")
                    call spec%disp%show(css_type([character(len(spec%domainAxisName%val),SKG) :: "quantile", spec%domainAxisName%val]), format = spec%reportFile%format%fixform, bmsize = 0_IK)
                    call spec%disp%skip(count = spec%disp%tmsize)
                    do iq = 1, size(stat%sample%quantile%prob, 1, IK)
                        write(spec%reportFile%unit, spec%reportFile%format%allreal) stat%sample%quantile%prob(iq), (stat%sample%quantile%quan(iq, idim), idim = 1, ndim)
                    end do
                    call spec%disp%skip(count = spec%disp%bmsize)
                    call spec%disp%note%show("This is the quantiles table of the variables of the final output refined sample.")

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Begin inter-chain convergence test in multiChain parallelization mode
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if                 CAF_ENABLED || MPI_ENABLED || OMP_ENABLED
                    blockInterChainConvergence: if (spec%parallelism%is%multiChain .and. spec%image%count > 1_IK) then

                        call spec%disp%note%show("Computing the inter-chain convergence probabilities...")
                        if (spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                            call spec%disp%note%show("Computing the inter-chain convergence probabilities...", unit = output_unit, bmsize = 0_IK)
                        end if
                        call spec%image%sync()

                        ! compute and report the KS convergence probabilities for all images.

                        multiChainConvergenceTest: block

                            real(RKG), allocatable :: sampleLogFuncState1(:,:)
                            real(RKG), allocatable :: sampleLogFuncState2(:,:)
                            integer(IK) :: imageID, idimMinProbKS, pidMinProbKS
                            character(:, SK), allocatable :: sampleFilePath
                            real(RKG), allocatable :: probKS(:)
                            real(RKG) :: minProbKS

                            minProbKS = 2._RKG !huge(minProbKS)
                            call setResized(probKS, ndim + 1_IK)
                            !call setRebound(probKS, 0_IK, ndim)

                            ! sort the refined chain on the current image.

                            err%stat = getErrTableRead(spec%sampleFile%file, sampleLogFuncState1, trans, sep = spec%outputSeparator%val, roff = 1_IK)
                            err%occurred = err%stat /= 0_IK
                            if (err%occurred) then
                                err%msg = PROCEDURE_NAME//SKG_": "//err%msg ! LCOV_EXCL_LINE
                                exit blockLeaderPostProcessing ! LCOV_EXCL_LINE
                            end if
                            do idim = 1, ndim + 1
                                call setSorted(sampleLogFuncState1(:, idim))
                            end do

                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            call spec%disp%show("stats.chain.refined.kstest.prob")
                            call spec%disp%show([css_type("processID"), stat%sfc%colname], format = spec%reportFile%format%fixform, bmsize = 0_IK)

                            do imageID = 1, spec%image%count

                                if (spec%image%id /= imageID) then

                                    ! read the refined chain on the other image.
                                    sampleFilePath = getReplaced(spec%sampleFile%file, spec%sampleFile%suffix, getReplaced(spec%sampleFile%suffix, SKG_"_pid"//getStr(spec%image%id), SKG_"_pid"//getStr(imageID)))
                                    err%stat = getErrTableRead(sampleFilePath, sampleLogFuncState2, trans, sep = spec%outputSeparator%val, roff = 1_IK)
                                    err%occurred = err%stat /= 0_IK
                                    if (err%occurred) then
                                        err%msg = PROCEDURE_NAME//SKG_": "//err%msg ! LCOV_EXCL_LINE
                                        exit blockLeaderPostProcessing ! LCOV_EXCL_LINE
                                    end if

                                    do idim = 1, ndim + 1
                                        ! sort the refined chain on the other image.
                                        call setSorted(sampleLogFuncState2(:, idim))
                                        ! compute the inter-chain KS probability table.
                                        probKS(idim) = getProbKS(getDisKolm(sampleLogFuncState1(:, idim), sampleLogFuncState2(:, idim), ascending), size(sampleLogFuncState1, 1, IK), size(sampleLogFuncState2, 1, IK))
                                        if (minProbKS <= probKS(idim)) cycle
                                        minProbKS = probKS(idim)
                                        pidMinProbKS = imageID
                                        idimMinProbKS = idim
                                    end do

                                    ! write the inter-chain KS probability table row

                                    write(spec%reportFile%unit, spec%reportFile%format%intreal) imageID, probKS

                                end if

                            end do

                            call spec%disp%skip(count = spec%disp%bmsize)
                            spec%msg = SKG_"This is the table of pairwise inter-chain Kolmogorov-Smirnov (KS) convergence (similarity) probabilities. &
                            &Higher KS probabilities are better, presenting less evidence for a lack of convergence of all chains to the same target density function."
                            call spec%disp%note%show(spec%msg)

                            ! write the smallest KS probabilities

                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            call spec%disp%show("stats.chain.refined.kstest.prob.min")
                            call spec%disp%show(minProbKS, format = format)
                            spec%msg = SKG_"This is the smallest KS-test probability for the inter-chain sampling convergence, which has happened between "//&
                            stat%sfc%colname(idimMinProbKS)%val//SKG_" on the chains generated by processes "//getStr(spec%image%id)//SKG_" and "//getStr(pidMinProbKS)//SKG_"."
                            call spec%disp%note%show(spec%msg)

                            ! Report the smallest KS probabilities on stdout.

                            if (.not. spec%outputSplashMode%is%silent) then
                                if (spec%image%is%first) then
                                    call spec%disp%note%show("The smallest KS probabilities for the inter-chain sampling convergence (higher is better):", unit = output_unit, bmsize = 0_IK)!, tmsize = 2_IK
                                end if
                                do imageID = 1, spec%image%count
                                    if (spec%image%id == imageID) then
                                        spec%msg = getStr(minProbKS)//SKG_" for "//stat%sfc%colname(idimMinProbKS)%val//&
                                        SKG_" on the chains generated by processes "//getStr(spec%image%id)//SKG_" and "//getStr(pidMinProbKS)//SKG_"."
                                        call spec%disp%note%show(spec%msg, unit = output_unit, tmsize = 0_IK, bmsize = 0_IK)
                                    end if
                                    flush(output_unit)
                                    call spec%image%sync()
                                end do
                                if (spec%image%is%first) then
                                    call execute_command_line("", cmdstat = err%stat)
                                    if (1 < spec%disp%bmsize) call spec%disp%skip(unit = output_unit, count = spec%disp%bmsize - 1)
                                end if
                            end if

                        end block multiChainConvergenceTest

                        call spec%image%sync()

                    end if blockInterChainConvergence
#endif
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! End inter-chain convergence test in multiChain parallelization mode
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end if blockSampleFileGeneration

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! End of generating the i.i.d. sample statistics and output file (if requested)
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! Mission accomplished.

                call spec%disp%note%show("Mission Accomplished.", tmsize = 2_IK * spec%disp%note%tmsize)!, tmsize = 3_IK, bmsize = 1_IK
                if (spec%reportFile%unit /= output_unit .and. spec%image%is%first .and. .not. spec%outputSplashMode%is%silent) then
                    flush(output_unit)
                    call execute_command_line("", cmdstat = err%stat)
                    call spec%disp%note%show("Mission Accomplished.", unit = output_unit, tmsize = 1_IK)!, tmsize = 1_IK, bmsize = 3_IK)
                end if
                close(spec%progressFile%unit)
                close(spec%reportFile%unit)
                err%msg = SK_""

            end associate

        end if blockLeaderPostProcessing

        ! A global sync is essential for parallel applications.

        SET_CAFMPI(call spec%image%sync())
        SET_CAFMPI(if (spec%parallelismMpiFinalizeEnabled%val) call spec%image%finalize())
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  RETURN_IF_FAILED
#undef  SET_DRAMDISE
#undef  SET_ParaNest
#undef  SET_CAFMPI
#undef  SET_SERIAL
#undef  SET_DEBUG
#undef  SET_CAF
#undef  SET_MPI
#undef  SET_OMP
#undef  COINDEX
#undef  TRACE