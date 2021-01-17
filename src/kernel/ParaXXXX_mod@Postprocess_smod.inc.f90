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
!> This file implements the body of the `Postprocess_smod` submodules of the ParaMonte `ParaDRAM_mod`, `ParaDISE_mod`, `ParaNest_mod` modules.
!>
!> \remark
!> This module requires preprocessing, prior to compilation.
!>
!> \author Amir Shahmoradi

#if !(defined PARADRAM || defined PARADISE || defined PARANEST)
#error "Unrecognized sampler in ParaXXXX_mod@Postprocess_smod.inc.f90"
#endif

!submodule (ParaDRAM_mod) Postprocess_smod
!submodule (ParaDISE_mod) Postprocess_smod
!submodule (ParaNest_mod) Postprocess_smod

    !use Constants_mod, only: IK, RK ! gfortran 9.3 compile crashes with this line
    implicit none

    character(*), parameter :: SUBMODULE_NAME = MODULE_NAME // "@Setup_smod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of 
    !> [ParaDRAM_type](@ref paradram_type)
    !> [ParaDISE_type](@ref paradise_type) 
    !> [ParaNest_type](@ref paranest_type) 
    !> classes.
    !> Perform the postprocessing analysis on the simulation results.
    !>
    !> @param[in]   self    :   An object instantiated from a ParaMonte sampler class.
    !>
    !> \remark
    !> This procedure requires preprocessing.
    module subroutine postprocess(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: postprocess
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit, stat_stopped_image
        use Decoration_mod, only: GENERIC_OUTPUT_FORMAT
        use Decoration_mod, only: GENERIC_TABBED_FORMAT
        use Decoration_mod, only: getGenericFormat
        use Decoration_mod, only: INDENT
        use Statistics_mod, only: getCorMatUpperFromCovMatUpper
        use Statistics_mod, only: getWeiSamCovUppMeanTrans
        use Statistics_mod, only: getQuantile
        use Statistics_mod, only: getMean
        use ParaMonte_mod, only: QPROB
        use Constants_mod, only: IK, RK, NLC, PMSM, UNDEFINED, POSINF_RK
        use DateTime_mod, only: getNiceDateTime
        use String_mod, only: num2str
        use Matrix_mod, only: getEye

        implicit none

#if defined PARADRAM
        class(ParaDRAM_type), intent(inout) :: self
#elif defined PARADISE
        class(ParaDISE_type), intent(inout) :: self
#elif defined PARANEST
        class(ParaNest_type), intent(inout) :: self
#endif
        character(*), parameter             :: PROCEDURE_NAME = SUBMODULE_NAME // "@runSampler()"

        integer(IK)                         :: i, iq, effectiveSampleSize
        character(:), allocatable           :: msg, formatStr, formatStrInt, formatStrReal, formatAllReal
        real(RK)    , allocatable           :: ContiguousChain(:,:) ! used to avoid temporary array creation and the compiler warning message in debug mode
        real(RK)                            :: mcmcSamplingEfficiency
        integer(IK)                         :: ndim

        ndim = self%nd%val

        ! @todo: Ever sinze postprocessing has become a separate submodule, is the following WARNING relevant anymore?
        ! WARNING: For any error that happens within the leader if-block below, the control must be passed to after the block
        ! WARNING: This is essential for error handling during testing in parallel mode.

        blockLeaderPostProcessing: if (self%Image%isLeader) then

            formatStrInt  = "('" // INDENT // "',1A" // self%LogFile%maxColWidth%str // ",*(I" // self%LogFile%maxColWidth%str // "))"
            formatStrReal = "('" // INDENT // "',1A" // self%LogFile%maxColWidth%str // ",*(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "E3))"

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.accepted"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%NumFunCall%accepted
            msg = "This is the total number of accepted function calls (unique samples)."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejected"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%NumFunCall%acceptedRejected
            msg = "This is the total number of accepted or rejected function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejectedDelayed"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%NumFunCall%acceptedRejectedDelayed
            msg = "This is the total number of accepted or rejected or delayed-rejection (if any requested) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.numFuncCall.acceptedRejectedDelayedUnused"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%NumFunCall%acceptedRejectedDelayedUnused
            msg = "This is the total number of accepted or rejected or unused function calls (by all processes, including delayed rejections, if any requested)."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.meanAcceptanceRate"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Chain%MeanAccRate(self%Stats%NumFunCall%accepted)
            msg = "This is the average MCMC acceptance rate."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            mcmcSamplingEfficiency = real(self%Stats%NumFunCall%accepted,kind=RK) / real(self%Stats%NumFunCall%acceptedRejected,kind=RK)

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejected"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) mcmcSamplingEfficiency
            msg =   "This is the MCMC sampling efficiency given the accepted and rejected function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayed"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(self%Stats%NumFunCall%accepted,kind=RK) / real(self%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
            msg =   "This is the MCMC sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayedUnused"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(self%Stats%NumFunCall%accepted,kind=RK) / real(self%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
            msg =   "This is the MCMC sampling efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls, that is, &
                    &the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.total"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Timer%Time%total
            msg = "This is the total runtime in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAccepted"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Timer%Time%total / real(self%Stats%NumFunCall%accepted,kind=RK)
            msg = "This is the average effective time cost of each accepted function call, in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejected"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Timer%Time%total / real(self%Stats%NumFunCall%acceptedRejected,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected function call, in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejectedDelayed"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Timer%Time%total / real(self%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected function call (including delayed-rejections, if any requested), in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCallAcceptedRejectedDelayedUnused"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Timer%Time%total / real(self%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
            msg = "This is the average effective time cost of each accepted or rejected or unused function call (including delayed-rejections, if any requested), in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (self%Image%count==1_IK) then
                msg = UNDEFINED ! LCOV_EXCL_LINE
            else ! LCOV_EXCL_LINE
                msg = num2str( self%Stats%avgCommTimePerFunCall ) ! LCOV_EXCL_LINE
            end if

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perInterProcessCommunication"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) msg
            msg = "This is the average time cost of inter-process communications per used (accepted or rejected or delayed-rejection) function call, in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.time.perFuncCall"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%avgTimePerFunCalInSec
            msg = "This is the average pure time cost of each function call, in seconds."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.current.numProcess"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Image%count
#if defined CAF_ENABLED || defined MPI_ENABLED
            if (self%Image%count==1_IK) then
                msg =   self%name // " " // &
                        "is being used in parallel mode but with only one processor. This is computationally inefficient. &
                        &Consider using the serial version of the code or provide more processes at runtime if it is beneficial."
                call self%note  ( prefix     = self%brand   &
                                , outputUnit = output_unit  &
                                , newline    = NLC          &
                                , marginTop  = 3_IK         &
                                , marginBot  = 0_IK         &
                                , msg        = msg          )
            end if
#endif
            msg = "This is the number of processes (images) used in this simulation."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin speedup compute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockSpeedup: block

                use Parallelism_mod, only: ForkJoin_type
                integer(IK)                 :: imageCount
                character(:), allocatable   :: formatIn, formatScaling
                logical                     :: isForkJoinParallelism
                type(ForkJoin_type)         :: ForkJoin
                character(:), allocatable   :: undefinedInfinity

                if (self%SpecBase%ParallelizationModel%isMultiChain) then
                    undefinedInfinity = "+INFINITY"
                else
                    undefinedInfinity = UNDEFINED
                end if

                isForkJoinParallelism = self%Image%count > 1_IK .and. self%SpecBase%ParallelizationModel%isSinglChain

                formatScaling = "('" // INDENT // "',10(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "E3))" ! ,:,','

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! compute the effective MCMC efficiency from the processor contributions and the current strong scaling
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                imageCount = 1_IK
                if (isForkJoinParallelism) imageCount = self%Image%count
                ForkJoin = ForkJoin_type( processCount = imageCount & ! LCOV_EXCL_LINE
                                        , lenProcessID = self%Stats%NumFunCall%accepted & ! LCOV_EXCL_LINE
                                        , ProcessID = self%Chain%ProcessID & ! LCOV_EXCL_LINE
                                        , successProb = mcmcSamplingEfficiency & ! LCOV_EXCL_LINE
                                        , seqSecTime = epsilon(1._RK) &  ! LCOV_EXCL_LINE ! time cost of the sequential section of the code, which is negligible here
                                        , parSecTime = self%Stats%avgTimePerFunCalInSec & ! LCOV_EXCL_LINE
                                        , comSecTime = self%Stats%avgCommTimePerFunCall & ! LCOV_EXCL_LINE
                                        )
                if (ForkJoin%Err%occurred) then
                    ! LCOV_EXCL_START
                    self%Err = ForkJoin%Err
                    self%Err%msg = PROCEDURE_NAME // self%Err%msg
                    call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                    exit blockLeaderPostProcessing
                    return
                    ! LCOV_EXCL_STOP
                end if

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.processContribution"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                block
                    integer, parameter :: NCOL = 40
                    character(:), allocatable :: formatInteger
                    formatInteger = "('"//INDENT//"',"//num2str(NCOL)//"(I0,:,' '))"
                    do imageCount = 1, ForkJoin%Contribution%count, NCOL
                        write(self%LogFile%unit,formatInteger) ForkJoin%Contribution%Frequency( imageCount : min( imageCount+NCOL-1, ForkJoin%Contribution%count ) )
                    end do
                end block
                msg =   "These are contributions of individual processes to the construction of the MCMC chain. &
                        &Essentially, they represent the total number of accepted states by the corresponding processor, &
                        &starting from the first processor to the last. This information is mostly informative in parallel &
                        &Fork-Join (singleChain) simulations."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.processContribution.geometricFit.successProbNormFac"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                if (isForkJoinParallelism) then
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%SuccessProb%PowellMinimum%xmin(1), exp(ForkJoin%SuccessProb%PowellMinimum%xmin(2))
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) UNDEFINED
                end if
                msg =   "These are the parameters of the Geometric fit to the distribution of the processor contributions &
                        &to the construction of the MCMC chain (the processor contributions are reported in the first column &
                        &of the output chain file. The fit has the following form: "//NLC//NLC// &
                        "    ProcessConstribution(i) = successProbNormFac(1) * successProbNormFac(2) * (1-successProbNormFac(1))^(i-1)"//NLC// &
                        "                            / (1 - (1 - successProbNormFac(1))^numProcess)"//NLC//NLC// &
                        "where i is the ID of the processor (starting from index 1), numProcess is the total number of processes &
                        &used in the simulation, and successProbNormFac(1) is equivalent to an effective MCMC sampling efficiency &
                        &computed from contributions of individual processes to the MCMC chain and successProbNormFac(2) is a &
                        &normalization constant."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.current.speedup"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%Speedup%current
                msg = "This is the estimated maximum speedup gained via "//self%SpecBase%ParallelizationModel%val//" parallelization model compared to serial mode."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.numProcess"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                if (isForkJoinParallelism) then
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%Speedup%Maximum%nproc
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) undefinedInfinity
                end if
                msg = "This is the predicted optimal number of physical computing processes for "//self%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.speedup"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                if (isForkJoinParallelism) then
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%Speedup%Maximum%value
                    msg = ""
                    if (ForkJoin%Speedup%current<1._RK) then
                        formatIn = "(g0.6)"
                        msg =   "The time cost of calling the user-provided objective function must be at least " // &
                                num2str(1._RK/ForkJoin%Speedup%current,formatIn) // " times more (that is, ~" // &
                                num2str(10**6*self%Stats%avgTimePerFunCalInSec/ForkJoin%Speedup%current,formatIn) // &
                                " microseconds) to see any performance benefits from " // &
                                self%SpecBase%ParallelizationModel%val // " parallelization model for this simulation. "
                        if (ForkJoin%Speedup%Maximum%nproc==1_IK) then
                            msg = msg// "Consider simulating in the serial mode in the future (if used within &
                                        &the same computing environment and with the same configuration as used here)."
                        else
                            msg = msg// "Consider simulating on " // num2str(ForkJoin%Speedup%Maximum%nproc) // " processors in the future &
                                        &(if used within the same computing environment and with the same configuration as used here)."
                        end if
                        call self%note  ( prefix   = self%brand     &
                                        , outputUnit = output_unit  &
                                        , newline    = NLC          &
                                        , marginTop  = 3_IK         &
                                        , marginBot  = 0_IK         &
                                        , msg        = msg          )
                        msg = NLC // msg
                    end if
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) undefinedInfinity
                end if
                msg = "This is the predicted optimal maximum speedup gained via "//self%SpecBase%ParallelizationModel%val//" parallelization model, given the current MCMC sampling efficiency." // msg
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.current.scaling.strong.speedup"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                if (isForkJoinParallelism) then
                    do imageCount = 1, ForkJoin%Speedup%count, 10
                        write(self%LogFile%unit,formatScaling) ForkJoin%Speedup%Scaling(imageCount:min(imageCount+9_IK,ForkJoin%Speedup%count))
                    end do
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) UNDEFINED
                end if
                msg =   "This is the predicted strong-scaling speedup behavior of the "//self%SpecBase%ParallelizationModel%val//" parallelization model, &
                        &given the current MCMC sampling efficiency, for increasing numbers of processes, starting from a single process."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! compute the absolute optimal parallelism efficiency under any MCMC sampling efficiency
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (isForkJoinParallelism) then
                    ForkJoin = ForkJoin_type( processCount = self%Image%count & ! LCOV_EXCL_LINE
                                            , lenProcessID = self%Stats%NumFunCall%accepted & ! LCOV_EXCL_LINE
                                            , ProcessID = self%Chain%ProcessID & ! LCOV_EXCL_LINE
                                            , successProb = 0._RK & ! LCOV_EXCL_LINE
                                            , seqSecTime = epsilon(1._RK) &  ! LCOV_EXCL_LINE ! time cost of the sequential section of the code, which is negligible here
                                            , parSecTime = self%Stats%avgTimePerFunCalInSec & ! LCOV_EXCL_LINE
                                            , comSecTime = self%Stats%avgCommTimePerFunCall & ! LCOV_EXCL_LINE
                                            )
                    if (ForkJoin%Err%occurred) then
                        ! LCOV_EXCL_START
                        self%Err = ForkJoin%Err
                        self%Err%msg = PROCEDURE_NAME // self%Err%msg
                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                        exit blockLeaderPostProcessing
                        return
                        ! LCOV_EXCL_STOP
                    end if
                end if ! otherwise, use the previously-generated ForkJoin

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.numProcess"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) ! LCOV_EXCL_LINE
                if (isForkJoinParallelism) then
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%Speedup%Maximum%nproc
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) undefinedInfinity
                end if
                msg = "This is the predicted absolute optimal number of physical computing processes for "//self%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.speedup"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) ! LCOV_EXCL_LINE
                if (isForkJoinParallelism) then
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) ForkJoin%Speedup%Maximum%value
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) undefinedInfinity
                end if
                msg =   "This is the predicted absolute optimal maximum speedup gained via "//self%SpecBase%ParallelizationModel%val//" parallelization model, under any MCMC sampling efficiency. &
                        &This simulation will likely NOT benefit from any additional computing processors beyond the predicted absolute optimal number, " // num2str(ForkJoin%Speedup%Maximum%nproc) // ", &
                        &in the above. This is true for any value of MCMC sampling efficiency. Keep in mind that the predicted absolute optimal number of processors is just an estimate &
                        &whose accuracy depends on many runtime factors, including the topology of the communication network being used, the number of processors per node, &
                        &and the number of tasks to each processor or node."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.parallelism.optimal.absolute.scaling.strong.speedup"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) ! LCOV_EXCL_LINE
                if (isForkJoinParallelism) then
                    do imageCount = 1, ForkJoin%Speedup%count, 10
                        write(self%LogFile%unit,formatScaling) ForkJoin%Speedup%Scaling(imageCount:min(imageCount+9_IK,ForkJoin%Speedup%count))
                    end do
                else
                    write(self%LogFile%unit,GENERIC_TABBED_FORMAT) UNDEFINED
                end if
                msg =   "This is the predicted absolute strong-scaling speedup behavior of the "//self%SpecBase%ParallelizationModel%val//" parallelization model, &
                        &under any MCMC sampling efficiency, for increasing numbers of processes, starting from a single process."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end block blockSpeedup


            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end speedup compute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.burnin.location.likelihoodBased"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%BurninLoc%compact
            msg = "This is the burnin location in the compact chain, based on the occurrence likelihood."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.burnin.location.adaptationBased"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%AdaptationBurninLoc%compact
            msg = "This is the burnin location in the compact chain, based on the value of burninAdaptationMeasure simulation specification."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.burnin.location.likelihoodBased"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%BurninLoc%verbose
            msg = "This is the burnin location in the verbose (Markov) chain, based on the occurrence likelihood."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.burnin.location.adaptationBased"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%AdaptationBurninLoc%verbose
            msg = "This is the burnin location in the verbose (Markov) chain, based on the value of burninAdaptationMeasure simulation specification."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! reset BurninLoc to the maximum value

            if (self%Stats%AdaptationBurninLoc%compact>self%Stats%BurninLoc%compact) then
                self%Stats%BurninLoc%compact = self%Stats%AdaptationBurninLoc%compact
                self%Stats%BurninLoc%verbose = self%Stats%AdaptationBurninLoc%verbose
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%LogFuncMode%val
            msg = "This is the maximum logFunc value (the maximum of the user-specified objective function)."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.compact.logFunc.max.location"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%LogFuncMode%Loc%compact
            msg = "This is the location of the first occurrence of the maximum logFunc in the compact chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max.location"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%LogFuncMode%Loc%verbose
            msg = "This is the location of the first occurrence of the maximum logFunc in the verbose (Markov) chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            formatAllReal = "('" // INDENT // "',*(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "E3))"

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.logFunc.max.state"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) (trim(adjustl(self%SpecBase%VariableNameList%Val(i))), i=1, ndim)
            write(self%LogFile%unit,formatAllReal) self%Stats%LogFuncMode%Crd
            msg = "This is the coordinates, within the domain of the user-specified objective function, where the maximum logFunc occurs."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the statistical properties of the MCMC chain
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (self%Image%isFirst) then
                call self%note  ( prefix        = self%brand    &
                                , outputUnit    = output_unit   &
                                , newline       = NLC           &
                                , marginTop     = 3_IK          &
                                , msg           = "Computing the statistical properties of the Markov chain..." )
            end if

            call self%Decor%writeDecoratedText  ( text = "\nThe statistical properties of the Markov chain\n" &
                                                , marginTop = 1_IK  &
                                                , marginBot = 1_IK  &
                                                , newline = "\n"    &
                                                , outputUnit = self%LogFile%unit )

            self%Stats%Chain%count = sum(self%Chain%Weight(self%Stats%BurninLoc%compact:self%Chain%count%compact))

            ! compute the covariance and correlation upper-triangle matrices

            if (allocated(self%Stats%Chain%Mean)) deallocate(self%Stats%Chain%Mean); allocate(self%Stats%Chain%Mean(ndim))
            if (allocated(self%Stats%Chain%CovMat)) deallocate(self%Stats%Chain%CovMat); allocate(self%Stats%Chain%CovMat(ndim,ndim))
            call getWeiSamCovUppMeanTrans   ( np = self%Chain%count%compact - self%Stats%BurninLoc%compact + 1_IK & ! LCOV_EXCL_LINE
                                            , sumWeight = self%Stats%Chain%count & ! LCOV_EXCL_LINE
                                            , nd = ndim & ! LCOV_EXCL_LINE
                                            , Point = self%Chain%State(1:ndim,self%Stats%BurninLoc%compact:self%Chain%count%compact) & ! LCOV_EXCL_LINE
                                            , Weight = self%Chain%Weight(self%Stats%BurninLoc%compact:self%Chain%count%compact) & ! LCOV_EXCL_LINE
                                            , CovMatUpper = self%Stats%Chain%CovMat & ! LCOV_EXCL_LINE
                                            , Mean = self%Stats%Chain%Mean & ! LCOV_EXCL_LINE
                                            )

            self%Stats%Chain%CorMat = getCorMatUpperFromCovMatUpper(nd=ndim,CovMatUpper=self%Stats%Chain%CovMat)

            ! transpose the covariance and correlation matrices

            do i = 1, ndim
                self%Stats%Chain%CorMat(i,i) = 1._RK
                self%Stats%Chain%CorMat(i+1:ndim,i) = self%Stats%Chain%CorMat(i,i+1:ndim)
                self%Stats%Chain%CovMat(i+1:ndim,i) = self%Stats%Chain%CovMat(i,i+1:ndim)
            end do

            ! compute the quantiles

            if (allocated(self%Stats%Chain%Quantile)) deallocate(self%Stats%Chain%Quantile)
            allocate(self%Stats%Chain%Quantile(QPROB%count,ndim))
            ContiguousChain = transpose(self%Chain%State(1:ndim,self%Stats%BurninLoc%compact:self%Chain%count%compact)) ! avoid temporary array creation and the warning message in debug mode
            do i = 1, ndim
                    self%Stats%Chain%Quantile(1:QPROB%count,i) = getQuantile( np = self%Chain%count%compact - self%Stats%BurninLoc%compact + 1_IK &
                                                                            , nq = QPROB%count & ! LCOV_EXCL_LINE
                                                                            , SortedQuantileProbability = QPROB%Value & ! LCOV_EXCL_LINE
                                                                            , Point = ContiguousChain(:,i) & ! LCOV_EXCL_LINE
                                                                            , Weight = self%Chain%Weight(self%Stats%BurninLoc%compact:self%Chain%count%compact) &
                                                                            , sumWeight = self%Stats%Chain%count & ! LCOV_EXCL_LINE
                                                                            )
            end do

            ! report the MCMC chain statistics

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.length.burninExcluded"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%Chain%count
            msg = "This is the length of the verbose (Markov) Chain excluding burnin."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.avgStd"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) "variableName", "Mean", "Standard Deviation"
            do i = 1, ndim
                write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Chain%Mean(i), sqrt(self%Stats%Chain%CovMat(i,i))
            end do
            msg = "This is the mean and standard deviation of the verbose (Markov) chain variables."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.covmat"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) "", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Chain%CovMat(1:ndim,i)
            end do
            msg = "This is the covariance matrix of the verbose (Markov) chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.cormat"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) "", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 1, ndim
                write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Chain%CorMat(1:ndim,i)
            end do
            msg = "This is the correlation matrix of the verbose (Markov) chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.verbose.quantile"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) "Quantile", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do iq = 1, QPROB%count
                write(self%LogFile%unit,formatStrReal) trim(adjustl(QPROB%Name(iq))), (self%Stats%Chain%Quantile(iq,i),i=1,ndim)
            end do
            msg = "This is the quantiles table of the variables of the verbose (Markov) chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Generate the i.i.d. sample statistics and output file (if requested)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! report refined sample statistics, and generate output refined sample if requested.

            if (self%Image%isFirst) then
                call self%note  ( prefix        = self%brand    &
                                , outputUnit    = output_unit   &
                                , newline       = NLC           &
                                , msg           = "Computing the final decorrelated sample size..." )
            end if

            call self%RefinedChain%get  ( CFC                       = self%Chain                                  &
                                        , Err                       = self%Err                                    &
                                        , burninLoc                 = self%Stats%BurninLoc%compact                &
                                        , sampleRefinementCount     = self%SpecMCMC%SampleRefinementCount%val     &
                                        , sampleRefinementMethod    = self%SpecMCMC%SampleRefinementMethod%val    &
                                        )

            if (self%Err%occurred) then
                ! LCOV_EXCL_START
                self%Err%msg = PROCEDURE_NAME // self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                exit blockLeaderPostProcessing
                return
                ! LCOV_EXCL_STOP
            end if

            ! compute the maximum integrated autocorrelation times for each variable

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            formatStr = "('" // INDENT // "',2I" // self%LogFile%maxColWidth%str // ",*(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "E3))"

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.iac"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,self%LogFile%format) "RefinementStage","SampleSize","IAC_SampleLogFunc",("IAC_"//trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
            do i = 0, self%RefinedChain%numRefinement
                write(self%LogFile%unit,formatStr) i, self%RefinedChain%Count(i)%Verbose, self%RefinedChain%IAC(0:ndim,i)
            end do
            msg = "This is the table of the Integrated Autocorrelation (IAC) of individual variables in the verbose (Markov) chain, at increasing stages of chain refinements."
            if (self%RefinedChain%numRefinement==0_IK) then
                msg = msg // NLC // "The user-specified sampleRefinementCount ("// num2str(self%SpecMCMC%SampleRefinementCount%val) // ") &
                                    &is too small to ensure an accurate computation of the decorrelated i.i.d. effective sample size. &
                                    &No refinement of the Markov chain was performed."
            end if
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! report the final Effective Sample Size (ESS) based on IAC

            !blockEffectiveSampleSize: associate( effectiveSampleSize => sum(self%RefinedChain%Weight(1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact)) )
            effectiveSampleSize = self%RefinedChain%Count(self%RefinedChain%numRefinement)%verbose

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.ess"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) effectiveSampleSize
            msg = "This is the estimated Effective (decorrelated) Sample Size (ESS) of the final refined chain."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAccepted"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(self%Stats%NumFunCall%accepted,kind=RK)
            msg =   "This is the effective MCMC sampling efficiency given the accepted function calls, that is, &
                    &the final refined effective sample size (ESS) divided by the number of accepted function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejected"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(self%Stats%NumFunCall%acceptedRejected,kind=RK)
            msg =   "This is the effective MCMC sampling efficiency given the accepted and rejected function calls, that is, &
                    &the final refined effective sample size (ESS) divided by the number of (accepted + rejected) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejectedDelayed"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(self%Stats%NumFunCall%acceptedRejectedDelayed,kind=RK)
            msg =   "This is the effective MCMC sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, that is, &
                    &the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.efficiency.essOverAcceptedRejectedDelayedUnused"
            write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
            write(self%LogFile%unit,GENERIC_TABBED_FORMAT) real(effectiveSampleSize,kind=RK) / real(self%Stats%NumFunCall%acceptedRejectedDelayedUnused,kind=RK)
            msg =   "This is the effective MCMC sampling efficiency given the accepted, rejected, delayed-rejection (if any requested), and unused function calls, that is, &
                    &the final refined effective sample size (ESS) divided by the number of (accepted + rejected + delayed-rejection + unused) function calls."
            call self%reportDesc(msg)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !end associate blockEffectiveSampleSize

            ! generate output refined sample if requested

            blockSampleFileGeneration: if (self%SpecBase%SampleSize%val==0_IK) then

                call self%note  ( prefix        = self%brand        & ! LCOV_EXCL_LINE
                                , outputUnit    = self%LogFile%unit & ! LCOV_EXCL_LINE
                                , newline       = NLC               & ! LCOV_EXCL_LINE
                                , msg           = "Skipping the generation of the decorrelated sample and output file, as requested by the user..." )

            else blockSampleFileGeneration

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! sample file generation report
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ! report to the report-file(s)

                call self%note  ( prefix        = self%brand        & ! LCOV_EXCL_LINE
                                , outputUnit    = self%LogFile%unit & ! LCOV_EXCL_LINE
                                , newline       = NLC               & ! LCOV_EXCL_LINE
                                , msg           = "Generating the output " // self%SampleFile%suffix // " file:"//NLC // self%SampleFile%Path%original )


                if (self%Image%isFirst) then

                    ! print the message for the generating the output sample file on the first image

                    call self%note  ( prefix     = self%brand       & ! LCOV_EXCL_LINE
                                    , outputUnit = output_unit      & ! LCOV_EXCL_LINE
                                    , newline    = NLC              & ! LCOV_EXCL_LINE
                                    , marginBot  = 0_IK             & ! LCOV_EXCL_LINE
                                    , msg        = "Generating the output " // self%SampleFile%suffix // " file:" )

                    call self%note  ( prefix     = self%brand       &
                                    , outputUnit = output_unit      &
                                    , newline    = NLC              &
                                    , marginTop  = 0_IK             &
                                    , marginBot  = 0_IK             &
                                    , msg = self%SampleFile%Path%original )

                    ! print the message for the generating the output sample file on the rest of the images in order

#if defined CAF_ENABLED || defined MPI_ENABLED
                    if (self%SpecBase%ParallelizationModel%isMultiChain) then
                        block
                            use String_mod, only: replaceStr, num2str
                            integer(IK) :: imageID
                            do imageID = 2, self%Image%count
                                call self%note  ( prefix = self%brand       &
                                                , outputUnit = output_unit  &
                                                , newline = NLC             &
                                                , marginTop = 0_IK          &
                                                , marginBot = 0_IK          &
                                                , msg = replaceStr( string = self%SampleFile%Path%original, search = "process_1", substitute = "process_"//num2str(imageID) ) )
                            end do
                        end block
                    end if
#endif
                    call self%Decor%write(output_unit,0,1)

                end if

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! begin sample file generation
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (self%SpecBase%SampleSize%val/=-1_IK) then

                    if (self%SpecBase%SampleSize%val<0_IK) then
                        self%SpecBase%SampleSize%abs = abs(self%SpecBase%SampleSize%abs) * self%RefinedChain%Count(self%RefinedChain%numRefinement)%verbose
                    end if

                    ! regenerate the refined sample, this time with the user-specified sample size.

                    call self%RefinedChain%get  ( CFC                       = self%Chain &
                                                , Err                       = self%Err &
                                                , burninLoc                 = self%Stats%BurninLoc%compact &
                                                , refinedChainSize          = self%SpecBase%SampleSize%abs &
                                                , sampleRefinementCount     = self%SpecMCMC%SampleRefinementCount%val     &
                                                , sampleRefinementMethod    = self%SpecMCMC%SampleRefinementMethod%val    &
                                                )
                    if (self%Err%occurred) then
                        ! LCOV_EXCL_START
                        self%Err%msg = PROCEDURE_NAME // self%Err%msg
                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                        exit blockLeaderPostProcessing
                        return
                        ! LCOV_EXCL_STOP
                    end if

                end if

                ! open the output sample file

                self%SampleFile%unit = 3000001 + self%Image%id  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
                open( unit      = self%SampleFile%unit            &
                    , file      = self%SampleFile%Path%original   &
                    , status    = self%SampleFile%status          &
                    , iostat    = self%SampleFile%Err%stat        &
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                    , SHARED &
#endif
                    , position  = self%SampleFile%Position%value  )
                self%Err = self%SampleFile%getOpenErr(self%SampleFile%Err%stat)
                if (self%Err%occurred) then
                    ! LCOV_EXCL_START
                    self%Err%msg = PROCEDURE_NAME//": Error occurred while opening the "//self%name//" "//self%SampleFile%suffix//" file='"//self%SampleFile%Path%original//"'. "
                    call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                    exit blockLeaderPostProcessing
                    return
                    ! LCOV_EXCL_STOP
                end if

                ! determine the sample file's contents' format

                if (self%SpecBase%OutputColumnWidth%val>0) then
                    formatStr = "(*(E"//self%SpecBase%OutputColumnWidth%str//"."//self%SpecBase%OutputRealPrecision%str//"E3,:,'"//self%SpecBase%OutputDelimiter%val//"'))"
                else
                    formatStr = self%SampleFile%format
                end if

                ! write to the output sample file

                call self%RefinedChain%write(self%SampleFile%unit,self%SampleFile%format,formatStr)
                close(self%SampleFile%unit)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Compute the statistical properties of the refined sample
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (self%Image%isFirst) then
                    call self%note  ( prefix        = self%brand        &
                                    , outputUnit    = output_unit       &
                                    , newline       = NLC               &
                                    , marginTop     = 2_IK              &
                                    , msg           = "Computing the statistical properties of the final refined sample..." )
                end if

                call self%Decor%writeDecoratedText  ( text = "\nThe statistical properties of the final refined sample \n" &
                                                    , marginTop = 1_IK  &
                                                    , marginBot = 1_IK  &
                                                    , newline = "\n"    &
                                                    , outputUnit = self%LogFile%unit )

                self%Stats%Sample%count = self%RefinedChain%Count(self%RefinedChain%numRefinement)%verbose

                ! compute the covariance and correlation upper-triangle matrices

                if (allocated(self%Stats%Sample%Mean)) deallocate(self%Stats%Sample%Mean); allocate(self%Stats%Sample%Mean(ndim))
                if (allocated(self%Stats%Sample%CovMat)) deallocate(self%Stats%Sample%CovMat); allocate(self%Stats%Sample%CovMat(ndim,ndim))
                if (self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact > ndim + 1_IK) then
                    ContiguousChain = self%RefinedChain%LogFuncState(1:ndim,1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact)
                    call getWeiSamCovUppMeanTrans   ( np = self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact &
                                                    , sumWeight = self%RefinedChain%Count(self%RefinedChain%numRefinement)%verbose &
                                                    , nd = ndim & ! LCOV_EXCL_LINE
                                                    , Point = ContiguousChain & ! LCOV_EXCL_LINE
                                                    , Weight = self%RefinedChain%Weight(1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact) & ! LCOV_EXCL_LINE
                                                    , CovMatUpper = self%Stats%Sample%CovMat & ! LCOV_EXCL_LINE
                                                    , Mean = self%Stats%Sample%Mean & ! LCOV_EXCL_LINE
                                                    )
                    self%Stats%Sample%CorMat = getCorMatUpperFromCovMatUpper(nd=ndim,CovMatUpper=self%Stats%Sample%CovMat)

                else
                    if (allocated(self%Stats%Sample%CovMat)) deallocate(self%Stats%Sample%CovMat) ! LCOV_EXCL_LINE
                    allocate(self%Stats%Sample%CovMat(ndim,ndim), source = POSINF_RK)
                    self%Stats%Sample%CorMat = getEye(ndim,ndim)
                    self%Stats%Sample%Mean   = getMean  ( nd = ndim &
                                                        , np = self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact &
                                                        , Point = self%RefinedChain%LogFuncState(1:ndim,1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact) & ! LCOV_EXCL_LINE
                                                        , Weight = self%RefinedChain%Weight(1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact) &
                                                        )
                end if

                ! transpose the covariance and correlation matrices

                do i = 1, ndim
                    self%Stats%Sample%CorMat(i,i) = 1._RK
                    self%Stats%Sample%CorMat(i+1:ndim,i) = self%Stats%Sample%CorMat(i,i+1:ndim)
                    self%Stats%Sample%CovMat(i+1:ndim,i) = self%Stats%Sample%CovMat(i,i+1:ndim)
                end do

                ! compute the quantiles

                ContiguousChain = transpose(self%RefinedChain%LogFuncState(1:ndim,1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact)) ! avoid temporary array creation and the warning message in debug mode
                if (allocated(self%Stats%Sample%Quantile)) deallocate(self%Stats%Sample%Quantile)
                allocate(self%Stats%Sample%Quantile(QPROB%count,ndim))
                do i = 1, ndim
                    self%Stats%Sample%Quantile(1:QPROB%count,i) = getQuantile   ( np = self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact &
                                                                                , nq = QPROB%count & ! LCOV_EXCL_LINE
                                                                                , SortedQuantileProbability = QPROB%Value & ! LCOV_EXCL_LINE
                                                                                , Point = ContiguousChain(:,i) & ! LCOV_EXCL_LINE
                                                                                , Weight = self%RefinedChain%Weight(1:self%RefinedChain%Count(self%RefinedChain%numRefinement)%compact) & ! LCOV_EXCL_LINE
                                                                                , sumWeight = self%RefinedChain%Count(self%RefinedChain%numRefinement)%verbose &
                                                                                )
                end do

                ! report the refined chain statistics

                !formatStr = "(1A" // self%LogFile%maxColWidth%str // ",*(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "))"

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.length"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_TABBED_FORMAT) self%Stats%Sample%count
                msg = "This is the final output refined sample size. "
                if (self%SpecBase%SampleSize%val/=-1_IK) then
                    if (abs(self%SpecBase%SampleSize%abs)<effectiveSampleSize) then
                        msg = msg // "The user-requested sample size ("// num2str(self%SpecBase%SampleSize%abs) // ") is smaller &
                            &than the potentially-optimal i.i.d. sample size (" // num2str(effectiveSampleSize) // "). The output sample &
                            &contains i.i.d. samples, however, the sample-size could have been larger if it had been set to the optimal size. &
                            &To get the optimal size in the future runs, set sampleSize = -1, or drop it from the input list."
                    elseif (abs(self%SpecBase%SampleSize%abs)>effectiveSampleSize) then
                        msg = msg // "The user-requested sample size ("// num2str(self%SpecBase%SampleSize%abs) // ") is larger than &
                            &the potentially-optimal i.i.d. sample size (" // num2str(effectiveSampleSize) // "). The resulting sample &
                            &likely contains duplicates and is not independently and identically distributed (i.i.d.).\nTo get the optimal &
                            &size in the future runs, set sampleSize = -1, or drop it from the input list."
                    else ! LCOV_EXCL_LINE
                        msg = msg // "How lucky that could be! The user-requested sample size (" // num2str(self%SpecBase%SampleSize%abs) // & ! LCOV_EXCL_LINE
                             ") is equal to the potentially-optimal i.i.d. sample size determined by the "//self%name//" sampler." ! LCOV_EXCL_LINE
                    end if
                end if
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.avgStd"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,self%LogFile%format) "variableName", "Mean", "Standard Deviation"
                do i = 1, ndim
                    write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Sample%Mean(i), sqrt(self%Stats%Sample%CovMat(i,i))
                end do
                msg = "This is the Mean and standard deviation table of the final output refined sample."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.covmat"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,self%LogFile%format) "", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Sample%CovMat(1:ndim,i)
                end do
                msg = "This is the covariance matrix of the final output refined sample."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.cormat"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,self%LogFile%format) "", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do i = 1, ndim
                    write(self%LogFile%unit,formatStrReal) trim(adjustl(self%SpecBase%VariableNameList%Val(i))), self%Stats%Sample%CorMat(1:ndim,i)
                end do
                msg = "This is the correlation matrix of the final output refined sample."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.quantile"
                write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                write(self%LogFile%unit,self%LogFile%format) "Quantile", (trim(adjustl(self%SpecBase%VariableNameList%Val(i))),i=1,ndim)
                do iq = 1, QPROB%count
                    write(self%LogFile%unit,formatStrReal) trim(adjustl(QPROB%Name(iq))), (self%Stats%Sample%Quantile(iq,i),i=1,ndim)
                end do
                msg = "This is the quantiles table of the variables of the final output refined sample."
                call self%reportDesc(msg)

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Begin inter-chain convergence test in multiChain parallelization mode
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined CAF_ENABLED || defined MPI_ENABLED

                blockInterChainConvergence: if (self%SpecBase%ParallelizationModel%isMultiChain .and. self%Image%count>1_IK) then

                    call self%note  ( prefix            = self%brand        &
                                    , outputUnit        = self%LogFile%unit &
                                    , newline           = NLC               &
                                    , marginTop         = 1_IK              &
                                    , marginBot         = 1_IK              &
                                    , msg               = "Computing the inter-chain convergence probabilities..." )
                    if (self%Image%isFirst) then ! print the message on stdout
                        call self%note  ( prefix        = self%brand        &
                                        , outputUnit    = output_unit       &
                                        , newline       = NLC               &
                                        , marginTop     = 1_IK              &
                                        , marginBot     = 1_IK              &
                                        , msg           = "Computing the inter-chain convergence probabilities..." )
                    end if

                    call self%Image%sync()

                    ! read the sample files generated by other images

                    multiChainConvergenceTest: block

                        use Sort_mod, only: sortAscending
                        use String_mod, only: replaceStr
                        use Statistics_mod, only: doSortedKS2
#if defined PARADRAM
                        use ParaDRAM_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
#elif defined PARADISE
                        use ParaDISE_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
#elif defined PARANEST
                        use ParaNest_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
#endif
                        type(RefinedChain_type)     :: RefinedChainThisImage, RefinedChainThatImage
                        integer(IK)                 :: imageID,indexMinProbKS,imageMinProbKS
                        real(RK)                    :: statKS, minProbKS
                        real(RK)    , allocatable   :: ProbKS(:)
                        character(:), allocatable   :: inputSamplePath

                        minProbKS = huge(minProbKS)
                        if (allocated(ProbKS)) deallocate(ProbKS) ! LCOV_EXCL_LINE
                        allocate(ProbKS(0:ndim))

                        ! read the refined chain on the current image

                        RefinedChainThisImage = readRefinedChain( sampleFilePath=self%SampleFile%Path%original, delimiter=self%SpecBase%OutputDelimiter%val, ndim=ndim, tenabled = .true. )
                        if (RefinedChainThisImage%Err%occurred) then
                            ! LCOV_EXCL_START
                            self%Err%occurred = .true.
                            self%Err%msg = PROCEDURE_NAME//RefinedChainThisImage%Err%msg
                            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                            exit blockLeaderPostProcessing
                            return
                            ! LCOV_EXCL_STOP
                        end if

                        ! sort the refined chain on the current image

                        do i = 0, ndim
                            call sortAscending  ( np = RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose &
                                                , Point = RefinedChainThisImage%LogFuncState(1:RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose,i) &
                                                , Err = self%Err &
                                                )
                            if (self%Err%occurred) then
                                ! LCOV_EXCL_START
                                self%Err%msg = PROCEDURE_NAME//self%Err%msg
                                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                                exit blockLeaderPostProcessing
                                return
                                ! LCOV_EXCL_STOP
                            end if
                        end do

                        ! compute and report the KS convergence probabilities for all images

                        formatStr = "(1I" // self%LogFile%maxColWidth%str // ",*(E" // self%LogFile%maxColWidth%str // "." // self%SpecBase%OutputRealPrecision%str // "E3))"

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.kstest.prob"
                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(self%LogFile%unit,self%LogFile%format) "ProcessID",("probKS("//trim(adjustl(self%RefinedChain%ColHeader(i)%record))//")",i=0,ndim)

                        do imageID = 1, self%Image%count

                            if (imageID/=self%Image%id) then

                                ! read the refined chain on the other image

                                inputSamplePath = replaceStr( string = self%SampleFile%Path%original &
                                                            , search = "process_"//num2str(self%Image%id) &
                                                            , substitute = "process_"//num2str(imageID) )

                                RefinedChainThatImage = readRefinedChain( sampleFilePath=inputSamplePath, delimiter=self%SpecBase%OutputDelimiter%val, ndim=ndim, tenabled = .true. )
                                if (RefinedChainThatImage%Err%occurred) then
                                    ! LCOV_EXCL_START
                                    self%Err%occurred = .true.
                                    self%Err%msg = PROCEDURE_NAME//RefinedChainThatImage%Err%msg
                                    call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                                    exit blockLeaderPostProcessing
                                    return
                                    ! LCOV_EXCL_STOP
                                end if

                                do i = 0, ndim

                                    ! sort the refined chain on the other image

                                    call sortAscending  ( np = RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose &
                                                        , Point = RefinedChainThatImage%LogFuncState(1:RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose,i) &
                                                        , Err = self%Err &
                                                        )
                                    if (self%Err%occurred) then
                                        ! LCOV_EXCL_START
                                        self%Err%msg = PROCEDURE_NAME//self%Err%msg
                                        call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                                        exit blockLeaderPostProcessing
                                        return
                                        ! LCOV_EXCL_STOP
                                    end if

                                    ! compute the inter-chain KS probability table

                                    call doSortedKS2( np1 = RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose &
                                                    , np2 = RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose &
                                                    , SortedPoint1 = RefinedChainThisImage%LogFuncState(1:RefinedChainThisImage%Count(RefinedChainThisImage%numRefinement)%verbose,i) &
                                                    , SortedPoint2 = RefinedChainThatImage%LogFuncState(1:RefinedChainThatImage%Count(RefinedChainThatImage%numRefinement)%verbose,i) &
                                                    , statKS = statKS &
                                                    , probKS = ProbKS(i) &
                                                    )

                                    ! determine the minimum KS

                                    if (ProbKS(i)<minProbKS) then
                                        minProbKS = ProbKS(i)
                                        indexMinProbKS = i
                                        imageMinProbKS = imageID
                                    end if

                                end do

                                ! write the inter-chain KS probability table row

                                write(self%LogFile%unit, formatStr) imageID, (ProbKS(i), i = 0, ndim)

                            end if

                        end do
                        msg =   "This is the table pairwise inter-chain Kolmogorov-Smirnov (KS) convergence (similarity) probabilities. &
                                &Higher KS probabilities are better, indicating less evidence for a lack of convergence."
                        call self%reportDesc(msg)

                        ! write the smallest KS probabilities

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT) "stats.chain.refined.kstest.prob.min"
                        write(self%LogFile%unit,GENERIC_OUTPUT_FORMAT)
                        write(self%LogFile%unit,GENERIC_TABBED_FORMAT) minProbKS
                        msg =   "This is the smallest KS-test probability for the inter-chain sampling convergence, which has happened between "//self%RefinedChain%ColHeader(indexMinProbKS)%record// &
                                " on the chains generated by processes "//num2str(self%Image%id)//" and "//num2str(imageMinProbKS)//"."
                        call self%reportDesc(msg)

                        ! report the smallest KS probabilities on stdout

                        if (self%Image%isFirst) call self%note  ( prefix     = self%brand       &
                                                                , outputUnit = output_unit      &
                                                                , newline    = NLC              &
                                                                , marginTop  = 2_IK             &
                                                                , marginBot  = 0_IK             &
                                                                , msg        = "The smallest KS probabilities for the inter-chain sampling convergence:" )

                        do imageID = 1, self%Image%count
                            if (imageID==self%Image%id) then
                                call self%note  ( prefix     = self%brand       &
                                                , outputUnit = output_unit      &
                                                , newline    = NLC              &
                                                , marginTop  = 0_IK             &
                                                , marginBot  = 0_IK             &
                                                , msg        = num2str(minProbKS)//" for "//self%RefinedChain%ColHeader(indexMinProbKS)%record//&
                                                               " on the chains generated by processes "//num2str(self%Image%id)//&
                                                               " and "//num2str(imageMinProbKS)//"." )
                            end if
#if defined CAF_ENABLED || MPI_ENABLED
                            call execute_command_line(" ", cmdstat = self%Err%stat)
                            flush(output_unit)
                            call self%Image%sync()
#endif
                        end do

                    end block multiChainConvergenceTest

                end if blockInterChainConvergence
#endif

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! End inter-chain convergence test in multiChain parallelization mode
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end if blockSampleFileGeneration

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! End of generating the i.i.d. sample statistics and output file (if requested)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call self%Decor%write(self%LogFile%unit,1,1)

            ! Mission accomplished.

            call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = "\n", msg = "Mission Accomplished." )
            if (output_unit/=self%LogFile%unit .and. self%Image%isFirst) then
                call self%Decor%write(output_unit,1,1)
                call self%note( prefix = self%brand, outputUnit = output_unit, newline = "\n", msg = "Mission Accomplished." )
                call self%Decor%write(output_unit,1,1)
            end if

            close(self%TimeFile%unit)
            close(self%LogFile%unit)

        end if blockLeaderPostProcessing

    end subroutine postprocess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!end submodule Postprocess_smod