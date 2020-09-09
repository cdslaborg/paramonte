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

module Parallelism_mod

    use Optimization_mod, only: PowellMinimum_type
    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@Parallelism_mod"

    type, private :: Maximum_type
        real(RK)                    :: value
        integer(IK)                 :: nproc
    end type Maximum_type

    type, private :: Speedup_type
        integer(IK)                 :: count
        real(RK)                    :: current
        type(Maximum_type)          :: Maximum
        real(RK)    , allocatable   :: Scaling(:)
    end type Speedup_type

    type, private :: UniqueProcess_type
        integer(IK)                 :: count
        integer(IK) , allocatable   :: Identity(:)
        integer(IK) , allocatable   :: Frequency(:)
    end type UniqueProcess_type

    type, private, extends(UniqueProcess_type) :: Contribution_type
        real(RK) , allocatable      :: LogFrequency(:)
    end type Contribution_type

    type, private :: SuccessProb_type
        real(RK)                    :: current
        real(RK)                    :: effective
        type(PowellMinimum_type)    :: PowellMinimum
    end type SuccessProb_type

    type :: ForkJoin_type
        type(UniqueProcess_type)    :: UniqueProcess
        type(Contribution_type)     :: Contribution ! similar to UniqueProcess, but including all processes, including zero contributors
        type(SuccessProb_type)      :: SuccessProb
        type(Speedup_type)          :: Speedup
        type(Err_type)              :: Err
    contains
        procedure   , nopass        :: getForkJoinSpeedup
    end type ForkJoin_type

    interface ForkJoin_type
        module procedure :: constructForkJoin
    end interface ForkJoin_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructForkJoin(processCount, lenProcessID, ProcessID, successProb, seqSecTime, parSecTime, comSecTime) result(ForkJoin) ! nonpure

        use Statistics_mod, only: fitGeoCyclicLogPDF
        use Constants_mod, only: RK, IK
        use String_mod, only: num2str
        use Misc_mod, only: findUnique
        use Sort_mod, only: indexArray

        implicit none

        integer(IK) , intent(in)    :: processCount, lenProcessID, ProcessID(lenProcessID)
        real(RK)    , intent(in)    :: successProb, seqSecTime, parSecTime, comSecTime
        real(RK)                    :: targetSuccessProb
        type(ForkJoin_type)         :: ForkJoin

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME // "@constructForkJoin()"
        real(RK)    , parameter     :: LOG_HALF = -0.693147180559945_RK

        integer(IK) , allocatable   :: Indx(:)
        integer(IK)                 :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! check for invalid input
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ForkJoin%Err%occurred = .false.

        if (successProb<0._RK) then
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME // ": successProb must be a number between zero and one. The input value is: " // num2str(successProb)
            return
        end if

        if (processCount==1_IK) then
            ForkJoin%Speedup%Maximum%value = 1._RK
            ForkJoin%Speedup%Maximum%nproc = 1_IK
            ForkJoin%Speedup%current = 1._RK
            ForkJoin%Speedup%count = 1_IK
            ForkJoin%Speedup%Scaling = [1._RK]
            ForkJoin%Contribution%count = processCount
            ForkJoin%Contribution%Identity = [1_IK]
            ForkJoin%Contribution%Frequency = [lenProcessID]
            ForkJoin%Contribution%LogFrequency = log(real(ForkJoin%Contribution%Frequency,kind=RK))
            ForkJoin%UniqueProcess%count = 1_IK
            ForkJoin%UniqueProcess%Identity = [1_IK]
            ForkJoin%UniqueProcess%Frequency = [lenProcessID]
            ForkJoin%SuccessProb%current = successProb
            ForkJoin%SuccessProb%effective = successProb
            return
        elseif (processCount<1_IK) then
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME // ": processCount cannot be less than one. The input value is: " // num2str(processCount)
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! compute the effective successProb (MCMC efficiency) from the Process contributions
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call findUnique ( lenVector = lenProcessID &
                        , Vector = ProcessID &
                        , UniqueValue = ForkJoin%UniqueProcess%Identity &
                        , UniqueCount = ForkJoin%UniqueProcess%Frequency &
                        , lenUnique = ForkJoin%UniqueProcess%count &
                        )
        !maxContributorProcessID = ForkJoin%UniqueProcess%Identity(maxloc(ForkJoin%UniqueProcess%Frequency))

        ! Sort ascending. This rearrangement may not be necessary anymore, since fitGeoCyclicLogPDF() is now order-agnostic.
        ! From the IO report perspective, however, it is important to have it ordered.

        allocate(Indx(ForkJoin%UniqueProcess%count))
        call indexArray( n = ForkJoin%UniqueProcess%count, Array = ForkJoin%UniqueProcess%Identity, Indx = Indx, Err = ForkJoin%Err )
        if (ForkJoin%Err%occurred) then
            ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
            return
        end if

        ! get all processes contributions

        ForkJoin%UniqueProcess%Identity(:) = ForkJoin%UniqueProcess%Identity(Indx)
        ForkJoin%UniqueProcess%Frequency(:) = ForkJoin%UniqueProcess%Frequency(Indx)
        deallocate(Indx)
        ForkJoin%Contribution%count = processCount
        if (allocated(ForkJoin%Contribution%LogFrequency)) deallocate(ForkJoin%Contribution%LogFrequency); allocate(ForkJoin%Contribution%LogFrequency(ForkJoin%Contribution%count), source = 0._RK)
        if (allocated(ForkJoin%Contribution%Frequency)) deallocate(ForkJoin%Contribution%Frequency); allocate(ForkJoin%Contribution%Frequency(ForkJoin%Contribution%count), source = 0_IK)
        if (allocated(ForkJoin%Contribution%Identity)) deallocate(ForkJoin%Contribution%Identity); allocate(ForkJoin%Contribution%Identity(ForkJoin%Contribution%count))
        ForkJoin%Contribution%Frequency(ForkJoin%UniqueProcess%Identity) = ForkJoin%UniqueProcess%Frequency
        ForkJoin%Contribution%Identity = [(i, i = 1, ForkJoin%Contribution%count)]
        if (ForkJoin%Contribution%Frequency(1)==0_IK) then
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME//": The contribution of the first process to the simulation is zero. This is highly unusual and requires further investigation."
            return
        end if
        do i = 1, ForkJoin%Contribution%count
            ! TODO: xxx a better solution instead of this ad-hoc approach would be to fit for the GeoCyclicCDF. Amir.
            if (ForkJoin%Contribution%Frequency(i)/=0_IK) then
                ForkJoin%Contribution%LogFrequency(i) = log(real(ForkJoin%Contribution%Frequency(i),kind=RK))
            else
                ForkJoin%Contribution%LogFrequency(i) = ForkJoin%Contribution%LogFrequency(i-1) + LOG_HALF
            end if
        end do

        ForkJoin%SuccessProb%current = successProb
        if (ForkJoin%SuccessProb%current>0._RK .and. ForkJoin%SuccessProb%current<1._RK) then

            ! fit for the effective successProb

            ForkJoin%SuccessProb%PowellMinimum = fitGeoCyclicLogPDF ( maxNumTrial = processCount &
                                                                    , numTrial = ForkJoin%Contribution%count &
                                                                    , SuccessStep = ForkJoin%Contribution%Identity &
                                                                    , LogCount = ForkJoin%Contribution%LogFrequency &
                                                                    )
            if (ForkJoin%SuccessProb%PowellMinimum%Err%occurred) then
                ForkJoin%Err = ForkJoin%SuccessProb%PowellMinimum%Err
                ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
                return
            end if

            targetSuccessProb = ForkJoin%SuccessProb%PowellMinimum%xmin(1)

        elseif (ForkJoin%SuccessProb%current==0._RK .or. ForkJoin%SuccessProb%current==1._RK) then

            targetSuccessProb = ForkJoin%SuccessProb%current

        else

            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME // ": The input successProb cannot be larger than one or less than zero. The current input value is "//num2str(ForkJoin%SuccessProb%current)
            return

        end if

        ! fit for the effective successProb

        call getForkJoinSpeedup ( successProb = targetSuccessProb &                             ! max(mcmcSamplingEfficiency, ForkJoin%SuccessProb%PowellMinimum%xmin(1)) & ! avoid unstable estimates of the effective efficiency.
                                , seqSecTime = epsilon(seqSecTime) &                            ! time cost of the sequential section of the code, which is negligible here.
                                , parSecTime = parSecTime &                                     ! the serial runtime for the parallel section of the code.
                                , comSecTimePerProc = comSecTime / (processCount - 1) &         ! the communication overhead for each additional image beyond master.
                                , minMaxNumProc = 2 * processCount &                            ! speedup will be computed at least up to this process count, if not more.
                                , Speedup = ForkJoin%Speedup%Scaling &                          ! returned Speedup.
                                , lenSpeedup = ForkJoin%Speedup%count &                         ! length of the returned Speedup vector.
                                , maxSpeedup = ForkJoin%Speedup%Maximum%value &                 ! returned maximum speedup.
                                , maxSpeedupNumProc = ForkJoin%Speedup%Maximum%nproc &          ! returned number of processes for maximum speedup.
                                , Err = ForkJoin%Err &
                                )
        if (ForkJoin%Err%occurred) then
            ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
            return
        end if
        ForkJoin%Speedup%current = ForkJoin%Speedup%Scaling(processCount)

    end function constructForkJoin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getForkJoinSpeedup(successProb, seqSecTime, parSecTime, comSecTimePerProc, minMaxNumProc, Speedup, lenSpeedup, maxSpeedupNumProc, maxSpeedup, Err)

        use Statistics_mod, only: getLogProbGeoCyclic
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        use Misc_mod, only: resize

        implicit none

        real(RK)        , intent(in)                :: successProb, seqSecTime, parSecTime, comSecTimePerProc
        integer(IK)     , intent(in)                :: minMaxNumProc ! must be at least 1
        real(RK)        , intent(out), allocatable  :: Speedup(:)
        integer(IK)     , intent(out)               :: lenSpeedup
        real(RK)        , intent(out)               :: maxSpeedup
        integer(IK)     , intent(out)               :: maxSpeedupNumProc
        type(Err_type)  , intent(out), optional     :: Err

        character(*)    , parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@constructForkJoin()"
        integer(IK)     , parameter                 :: ABS_MAX_NUM_PROC = 1000000_IK

        real(RK)                                    :: FirstImageContribution(1)
        real(RK)                                    :: serialRuntime, parSecTimePerProc, comSecTime
        integer(IK)                                 :: numProc, lenSpeedupNew !, maxNumProc
        logical                                     :: maxSpeedupFound, isPresentErr

        isPresentErr = present(Err)
        if (isPresentErr) Err%occurred = .false.

        lenSpeedup = minMaxNumProc
        allocate(Speedup(lenSpeedup))

        numProc = 1
        Speedup(numProc) = 1._RK
        maxSpeedup = Speedup(numProc)
        maxSpeedupNumProc = numProc
        serialRuntime = seqSecTime + parSecTime ! serial runtime of the code per function call

        maxSpeedupFound = .false.
        loopOverProc: do

            numProc = numProc + 1

            ! resize allocation if needed

            if (numProc>lenSpeedup) then
                if (maxSpeedupFound) exit loopOverProc
                lenSpeedupNew = 2 * lenSpeedup
                call resize(Speedup, from = lenSpeedup, to = lenSpeedupNew)
                lenSpeedup = lenSpeedupNew
            end if


            ! compute the fraction of work done by the first image

            if (successProb==0._RK) then
                FirstImageContribution(1) = 1._RK / numProc
            else
                FirstImageContribution = exp(getLogProbGeoCyclic( successProb = successProb &
                                                                , maxNumTrial = numProc &
                                                                , numTrial = 1_IK &
                                                                , SuccessStep = [1_IK] &
                                                                ) )
            end if

            ! effective runtime of the parallel-section of the code, when executed in parallel on numProc processes

            parSecTimePerProc = parSecTime * FirstImageContribution(1)

            ! communication time grows linearly with the number of processes.

            comSecTime = (numProc-1_IK) * comSecTimePerProc

            Speedup(numProc) = serialRuntime / (seqSecTime + parSecTimePerProc + comSecTime)
            if (maxSpeedup < Speedup(numProc)) then
                maxSpeedup = Speedup(numProc)
                maxSpeedupNumProc = numProc
            else
                maxSpeedupFound = .true.
            end if

            if (numProc<ABS_MAX_NUM_PROC) cycle loopOverProc

            if (isPresentErr) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // &
                            ": Failed to find the number of processes with which the maximum speedup occurs. The search continued up to " // &
                            num2str(ABS_MAX_NUM_PROC) // " processes."
            end if
            return

        end do loopOverProc

    end subroutine getForkJoinSpeedup

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Parallelism_mod