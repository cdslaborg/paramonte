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

!>  \brief This module contains procedures for computing the parallel performance of the parallel algorithms.
!>  \author Amir Shahmoradi

module Parallelism_mod

    use Optimization_mod, only: PowellMinimum_type
    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@Parallelism_mod"

    !> The [Image_type](@ref image_type) class.
    type                            :: Image_type
        integer(IK)                 :: id                   !<  The ID of the runtime parallel process: 1, 2, 3, ...
        integer(IK)                 :: count                !<  The total number of runtime parallel processes available.
        logical                     :: isFirst = .false.    !<  A logical flag indicating whether the current process is ID #1.
        logical                     :: isNotFirst = .false. !<  A logical flag indicating whether the current process is NOT ID #1.
        logical                     :: isLeader = .false.   !<  A logical flag indicating whether the current process is a leader. This must be defined by the user.
        logical                     :: isRooter = .false.   !<  A logical flag indicating whether the current process is a follower. This must be defined by the user.
        character(:), allocatable   :: name                 !<  The name of the current process as a string with the format `"@process(id)"`.
    contains
        procedure, nopass           :: sync => syncImages
        procedure, pass             :: query => queryImage
        procedure, nopass           :: finalize => finalizeImages
    end type Image_type

    type, private :: Maximum_type
        real(RK)                    :: value            !<  The maximum speedup attained (or attainable).
        integer(IK)                 :: nproc            !<  The required number of processes for maximum speedup.
    end type Maximum_type

    type, private :: Speedup_type
        integer(IK)                 :: count            !<  The size of the `Scaling` vector.
        real(RK)                    :: current          !<  The speedup given the current available number of processes.
        type(Maximum_type)          :: Maximum          !<  An object of type [Maximum_type](@ref maximum_type), containing the predicted maximum speedup and process count.
        real(RK)    , allocatable   :: Scaling(:)       !<  A real vector of length `Speedup_type%count` containing the predicted speedup for a range process counts.
    end type Speedup_type

    type, private :: UniqueProcess_type
        integer(IK)                 :: count            !<  The sizes of the two vector components `Identity` and `Frequency` of `UniqueProcess_type`.
        integer(IK) , allocatable   :: Identity(:)      !<  A vector of size `UniqueProcess_type%count` containing the unique IDs of processes, i.e., the ranks of processes, starting from 1.
        integer(IK) , allocatable   :: Frequency(:)     !<  The frequency with which the process IDs have contributed to the simulation at hand.
    end type UniqueProcess_type

    type, private, extends(UniqueProcess_type) :: Contribution_type
        real(RK) , allocatable      :: LogFrequency(:)  !<  The natural logarithm of the `Frequency` vector component of the type [UniqueProcess_type](@ref uniqueprocess_type).
    end type Contribution_type

    type, private :: SuccessProb_type
        real(RK)                    :: current          !<  The probability of success in the Bernoulli trial. For example, the sampling efficiency of the ParaDRAM sampler.
        real(RK)                    :: effective        !<  The computed effective probability of success inferred from fitting the Cyclic Geometric distribution.
        type(PowellMinimum_type)    :: PowellMinimum    !<  An object of class [PowellMinimum_type](@ref optimization_mod::powellminimum_type) containing
                                                        !<  information about the Cyclic Geometric fit to process contribution data.
    end type SuccessProb_type

    !> The `ForkJoin_type` class.
    type :: ForkJoin_type
        type(UniqueProcess_type)    :: UniqueProcess
        type(Contribution_type)     :: Contribution     !<  Contains information similar to `UniqueProcess`, but including all processes, even zero contributors.
        type(SuccessProb_type)      :: SuccessProb
        type(Speedup_type)          :: Speedup
        type(Err_type)              :: Err
    contains
        procedure   , nopass        :: getForkJoinSpeedup
    end type ForkJoin_type

    interface ForkJoin_type
        module procedure :: constructForkJoin
    end interface ForkJoin_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the statistics of the parallel processors available, depending on the type of parallelism requested.
    !> This is a dynamic member of the [Image_type](@ref image_type) class.
    !>
    !> @param[inout]    Image       :   An object of class [Image_type](@ref image_type).
    !>                                  On output, all properties of `Image` will be reset,
    !>                                  except the attributes `isLeader` and `isRooter`.
    !>
    !> \warning
    !> This routine must not contain any synchronization statements under any circumstances.
    subroutine queryImage(Image)

        use Constants_mod, only: RK, IK
        use String_mod, only: num2str
        implicit none
        class(Image_type), intent(inout) :: Image

        ! setup general processor / coarray image variables

#if defined CAF_ENABLED
        Image%id             = this_image()
        Image%count          = num_images()
#elif defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isInitialized
            call mpi_initialized( isInitialized, ierrMPI )
            if (.not. isInitialized) call mpi_init(ierrMPI)
            call mpi_comm_rank(mpi_comm_world, Image%id, ierrMPI)
            call mpi_comm_size(mpi_comm_world, Image%count, ierrMPI)
            Image%id = Image%id + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
        end block
#else
        Image%id            = 1_IK
        Image%count         = 1_IK
#endif
        Image%name          = "@process(" // num2str(Image%id) // ")"
        Image%isFirst       = Image%id==1_IK
        Image%isNotFirst    = Image%id/=1_IK
        !Image%isLeader      = .false.   ! ATTN: this is to be set by the user at runtime, depending on the type of parallelism.
        !Image%isRooter      = .false.   ! ATTN: this is to be set by the user at runtime, depending on the type of parallelism.

    end subroutine queryImage

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Synchronize all existing parallel images and return nothing.
    !> This is a static member of the [Image_type](@ref image_type) class.
    !>
    !> \warning
    !> This routine global Coarray and MPI synchronization barriers and
    !> therefore, must be called by all processes in the current simulation.
    subroutine syncImages()
        implicit none
#if defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            !logical     :: isFinalized
            !call mpi_finalized( isFinalized, ierrMPI )
            !if (.not. isFinalized) then
                call mpi_barrier(mpi_comm_world,ierrMPI)
            !end if
        end block
#elif defined CAF_ENABLED
        sync all
#endif
    end subroutine syncImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Finalize the current parallel simulation and return nothing.
    !> This is a static member of the [Image_type](@ref image_type) class.
    !>
    !> \warning
    !> MPI communications will be shut down upon calling this routine and further interprocess communications will be impossible.
    subroutine finalizeImages() ! LCOV_EXCL_LINE
#if defined CAF_ENABLED
        sync all
#elif defined MPI_ENABLED
        use mpi
        implicit none
        integer(IK) :: ierrMPI
        logical     :: isFinalized
        call mpi_finalized( isFinalized, ierrMPI )
        if (.not. isFinalized) then
            call mpi_barrier(mpi_comm_world,ierrMPI)
            call mpi_finalize(ierrMPI)
        end if
#endif
    end subroutine finalizeImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is the constructor of the [ForkJoin_type](@ref forkjoin_type) class.
    !> Return the predicted speedup of the parallel simulation given the input characteristics and timing information of the simulation.
    !>
    !> @param[in]   processCount    :   The number of processes in the simulation.
    !> @param[in]   lenProcessID    :   The length of `ProcessID` vector.
    !> @param[in]   ProcessID       :   The vector of process IDs.
    !> @param[in]   successProb     :   The success probability (the effective acceptance rate per objective function call).
    !> @param[in]   seqSecTime      :   The timing of the sequential sections of the code.
    !> @param[in]   parSecTime      :   The timing of the parallel sections of the code.
    !> @param[in]   comSecTime      :   The timing of the communication sections of the code.
    !>
    !> \return
    !> `ForkJoin` : An object of class [ForkJoin_type](@ref forkjoin_type) containing the parallelization speedup.
    function constructForkJoin(processCount, lenProcessID, ProcessID, successProb, seqSecTime, parSecTime, comSecTime) result(ForkJoin) ! nonpure

        use GeoCyclicFit_mod, only: fitGeoCyclicLogPDF ! LCOV_EXCL_LINE
        use Constants_mod, only: IK, RK, SQRT_EPS_RK
        use String_mod, only: num2str
        use Misc_mod, only: findUnique
        use Sort_mod, only: indexArray

        implicit none

        integer(IK) , intent(in)    :: processCount, lenProcessID, ProcessID(lenProcessID)
        real(RK)    , intent(in)    :: successProb, seqSecTime, parSecTime, comSecTime
        type(ForkJoin_type)         :: ForkJoin

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME // "@constructForkJoin()"
        real(RK)    , parameter     :: LOG_HALF = -0.693147180559945_RK

        integer(IK) , allocatable   :: Indx(:)
        integer(IK)                 :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! check for invalid input
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ForkJoin%Err%occurred = .false.

        if (successProb<-SQRT_EPS_RK .or. successProb>1._RK+SQRT_EPS_RK) then
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME // ": successProb must be a number between zero and one. The input value is: " // num2str(successProb)
            return
        end if

        if (processCount<1_IK) then
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME // ": processCount cannot be less than one. The input value is: " // num2str(processCount)
            return
        end if

        if (processCount==1_IK .or. successProb <= SQRT_EPS_RK .or. successProb >= 1._RK-SQRT_EPS_RK) then
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
            ! LCOV_EXCL_START
            ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
            return
            ! LCOV_EXCL_STOP
        end if

        ! get all processes contributions

        ForkJoin%UniqueProcess%Identity(:) = ForkJoin%UniqueProcess%Identity(Indx)
        ForkJoin%UniqueProcess%Frequency(:) = ForkJoin%UniqueProcess%Frequency(Indx)

        deallocate(Indx)

        ForkJoin%Contribution%count = processCount
        if (allocated(ForkJoin%Contribution%LogFrequency)) deallocate(ForkJoin%Contribution%LogFrequency);
        allocate(ForkJoin%Contribution%LogFrequency(ForkJoin%Contribution%count), source = 0._RK)
        if (allocated(ForkJoin%Contribution%Frequency)) deallocate(ForkJoin%Contribution%Frequency);
        allocate(ForkJoin%Contribution%Frequency(ForkJoin%Contribution%count), source = 0_IK)
        if (allocated(ForkJoin%Contribution%Identity)) deallocate(ForkJoin%Contribution%Identity);
        allocate(ForkJoin%Contribution%Identity(ForkJoin%Contribution%count))

        ForkJoin%Contribution%Frequency(ForkJoin%UniqueProcess%Identity) = ForkJoin%UniqueProcess%Frequency
        ForkJoin%Contribution%Identity = [(i, i = 1, ForkJoin%Contribution%count)]
        if (ForkJoin%Contribution%Frequency(1)==0_IK) then
            ! LCOV_EXCL_START
            ForkJoin%Err%occurred = .true.
            ForkJoin%Err%msg = PROCEDURE_NAME//": The contribution of the first process to the simulation is zero. This is highly unusual and requires further investigation."
            return
            ! LCOV_EXCL_STOP
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

        ! fit for the effective successProb

        ForkJoin%SuccessProb%PowellMinimum = fitGeoCyclicLogPDF ( maxNumTrial = processCount & ! LCOV_EXCL_LINE
                                                                , numTrial = ForkJoin%Contribution%count & ! LCOV_EXCL_LINE
                                                                , SuccessStep = ForkJoin%Contribution%Identity & ! LCOV_EXCL_LINE
                                                                , LogCount = ForkJoin%Contribution%LogFrequency & ! LCOV_EXCL_LINE
                                                                )
        if (ForkJoin%SuccessProb%PowellMinimum%Err%occurred) then
            ! LCOV_EXCL_START
            ForkJoin%Err = ForkJoin%SuccessProb%PowellMinimum%Err
            ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
            return
            ! LCOV_EXCL_STOP
        end if

        ForkJoin%SuccessProb%effective = ForkJoin%SuccessProb%PowellMinimum%xmin(1)

        ! fit for the effective successProb

        call getForkJoinSpeedup ( successProb = ForkJoin%SuccessProb%effective &                ! max(mcmcSamplingEfficiency, ForkJoin%SuccessProb%PowellMinimum%xmin(1)) & ! avoid unstable estimates of the effective efficiency.
                                , seqSecTime = epsilon(seqSecTime) &                            ! time cost of the sequential section of the code, which is negligible here.
                                , parSecTime = parSecTime &                                     ! the serial runtime for the parallel section of the code.
                                , comSecTimePerProc = comSecTime / (processCount - 1) &         ! the communication overhead for each additional image beyond master.
                                , minMaxNumProc = 2 * processCount &                            ! speedup will be computed at least up to this process count, if not more.
                                , Speedup = ForkJoin%Speedup%Scaling &                          ! returned Speedup.
                                , lenSpeedup = ForkJoin%Speedup%count &                         ! length of the returned `Speedup%Scaling` vector.
                                , maxSpeedup = ForkJoin%Speedup%Maximum%value &                 ! returned maximum speedup.
                                , maxSpeedupNumProc = ForkJoin%Speedup%Maximum%nproc &          ! returned number of processes for maximum speedup.
                                , Err = ForkJoin%Err &
                                )
        if (ForkJoin%Err%occurred) then
            ! LCOV_EXCL_START
            ForkJoin%Err%msg = PROCEDURE_NAME // ForkJoin%Err%msg
            return
            ! LCOV_EXCL_STOP
        end if
        ForkJoin%Speedup%current = ForkJoin%Speedup%Scaling(processCount)

    end function constructForkJoin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Predict the parallel simulation speedup for a range of possible processor counts.
    !>
    !> @param[in]   successProb         :   The success probability (the effective acceptance rate per objective function call).
    !> @param[in]   seqSecTime          :   The timing of the sequential sections of the code.
    !> @param[in]   parSecTime          :   The timing of the parallel sections of the code.
    !> @param[in]   comSecTimePerProc   :   The timing of the communication sections of the code per processor.
    !> @param[in]   minMaxNumProc       :   The minimum number of processes for which the speedup will be computed. It must be at least 1.
    !> @param[out]  Speedup             :   The vector of speedup values for different counts of processes.
    !> @param[out]  lenSpeedup          :   The length of the vector `Speedup`.
    !> @param[out]  maxSpeedupNumProc   :   The number of processes at which the maximum speedup happens.
    !> @param[out]  maxSpeedup          :   The maximum speedup for any number of processors.
    !> @param[out]  Err                 :   An object of [Err_type](@ref err_mod::err_type) class indicating whether and error has occurred upon return.
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
        integer(IK), parameter                      :: SuccessStep(*) = [1_IK]

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
                                                                , SuccessStep = SuccessStep &
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

            if (numProc<ABS_MAX_NUM_PROC) cycle loopOverProc ! LCOV_EXCL_LINE

            ! LCOV_EXCL_START
            if (isPresentErr) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // &
                            ": Failed to find the number of processes with which the maximum speedup occurs. The search continued up to " // &
                            num2str(ABS_MAX_NUM_PROC) // " processes."
            end if
            return
            ! LCOV_EXCL_STOP

        end do loopOverProc

    end subroutine getForkJoinSpeedup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Parallelism_mod ! LCOV_EXCL_LINE