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

module RandomSeed_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter :: MODULE_NAME = "@RandomSeed_mod"

    public
    private :: setRandomSeed, getRandomSeed

    type :: RandomSeed_type
        integer(IK)               :: size = -huge(1_IK)
        integer(IK)               :: imageID = -huge(1_IK)
        integer(IK), allocatable  :: Value(:)
        logical                   :: isRepeatable = .false.
        logical                   :: isImageDistinct = .true.
        character(:), allocatable :: info
        type(Err_type)            :: Err
    contains
        procedure, public :: set => setRandomSeed
        procedure, public :: get => getRandomSeed
    end type RandomSeed_type

    interface RandomSeed_type
        module procedure :: constructRandomSeed
    end interface RandomSeed_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    function constructRandomSeed(imageID, inputSeed, isRepeatable, isImageDistinct) result(RandomSeed)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructRandomSeed
#endif
        implicit none
        integer(IK) , intent(in)            :: imageID
        integer(IK) , intent(in), optional  :: inputSeed
        logical     , intent(in), optional  :: isRepeatable, isImageDistinct
        type(RandomSeed_type)               :: RandomSeed

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@constructRandomSeed()"

        RandomSeed%Err%occurred = .false.
        RandomSeed%Err%msg = ""
        RandomSeed%info = ""

        RandomSeed%imageID = imageID
        if (RandomSeed%imageID<1_IK) then
            RandomSeed%Err%occurred = .true.
            RandomSeed%Err%msg = PROCEDURE_NAME // ": Internal error occurred. imageID cannot be less than 1."
            return
        end if

        RandomSeed%isRepeatable = .false.
        if (present(isRepeatable)) RandomSeed%isRepeatable = isRepeatable

        RandomSeed%isImageDistinct = .true.
        if (present(isImageDistinct)) RandomSeed%isImageDistinct = isImageDistinct

        call RandomSeed%set(inputSeed)
        if (RandomSeed%Err%occurred) then
            RandomSeed%Err%msg = PROCEDURE_NAME // RandomSeed%Err%msg
            return
        end if

        call RandomSeed%get()

    end function constructRandomSeed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    subroutine getRandomSeed(RandomSeed)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandomSeed
#endif
        implicit none
        class(RandomSeed_type), intent(inout)   :: RandomSeed
        RandomSeed%Err%occurred = .false.
        RandomSeed%Err%msg = ""
        !if (allocated(RandomSeed%Value)) deallocate(RandomSeed%Value)
        if (.not. allocated(RandomSeed%Value)) then
            call random_seed(size = RandomSeed%size)
            allocate(RandomSeed%Value(RandomSeed%size))
        end if
        call random_seed(get = RandomSeed%Value)
    end subroutine getRandomSeed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    subroutine setRandomSeed(RandomSeed,inputSeed)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setRandomSeed
#endif
        use Constants_mod, only: IK, RK, HUGE_IK
        use iso_fortran_env, only: int64
        implicit none
        class(RandomSeed_type), intent(inout)   :: RandomSeed
        integer(IK), intent(in), optional       :: inputSeed
        integer(IK)                             :: offsetImageRandomSeed, i, scalarSeed
        integer(IK)                             :: values(8)

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@setRandomSeed()"

        RandomSeed%Err%occurred = .false.
        RandomSeed%Err%msg = ""

        call random_seed(size = RandomSeed%size)
        if ( allocated(RandomSeed%Value) ) deallocate(RandomSeed%Value)
        allocate( RandomSeed%Value(RandomSeed%size) )

        if (present(inputSeed)) then
            scalarSeed = abs(inputSeed)
        elseif (RandomSeed%isRepeatable) then
            scalarSeed = 12357913_IK   ! set the seed to something fixed so that all random number sequences can be regenerated
        else    ! simulation is not repeatable, initialize the seed to something random, different on each images
            call date_and_time(values=values)
            scalarSeed = abs(sum(values))
            do
                if (scalarSeed<=huge(scalarSeed) ) exit
                scalarSeed = scalarSeed - huge(scalarSeed)
            end do
            if (scalarSeed==0) then
                RandomSeed%Err%occurred = .true.
                RandomSeed%Err%msg = PROCEDURE_NAME // ": Random seed cannot be zero."
                return
            end if
        end if

        ! now use scalarSeed to construct the random seed on all images
        if (RandomSeed%isImageDistinct) then
            offsetImageRandomSeed = 127_IK * RandomSeed%size * (RandomSeed%imageID-1)
        else
            offsetImageRandomSeed = 0
        end if
        do i = 1, RandomSeed%size
            RandomSeed%Value(i) = HUGE_IK - scalarSeed - offsetImageRandomSeed - 127_IK * (i-1)
            if (RandomSeed%Value(i)<0_IK) then
                RandomSeed%Value(i) = -RandomSeed%Value(i) 
            else
                RandomSeed%Value(i) = HUGE_IK - RandomSeed%Value(i) 
            end if 
        end do
        call random_seed(put=RandomSeed%Value)

!block
!write(*,"(*(g0,:,' '))") 
!write(*,"(*(g0,:,' '))") "RandomSeed%Value", RandomSeed%Value
!write(*,"(*(g0,:,' '))") 
!end block


        ! ATTN: xxx Intel compilers - for some unknown reason, the first generated random number seems to be garbage
        ! so here, the random number generator is iterated a couple of times before further usage.
        ! This needs to be taken care of, in the future. This problem showed itself when StartPoint in ParaDRAM sampler were to be set randomly. 
        ! This is where the first instance of random number usage occurs in ParaDRAM sampler.
        ! write(*,*) "RandomSeedObj%imageID, co_RandomSeed(1)%Value(:): ", RandomSeedObj%imageID, co_RandomSeed(1)%Value(:)

        ! ATTN: A follow-up on the above issue with the Intel compiler which seems to be a compiler bug: In a truly bizarre behavior,
        ! the Intel compiler random numbers as generated by call random_number() in the Statistics_mod module, for example when called from
        ! ParaDRAMProposal_mod.inc.f90, are not repeatable even after reseting the random_seed. Even more bizarre is the observation that the 
        ! repeatability of the random numbers depends on the loop length (for example as implemented in the debugging of getRandGaus().
        ! The same behavior is also observed below, where any loop length less than ~30 yields non-repeatable random number sequences.
        ! This needs an in-depth investigation. Update: Such behavior was also observed with the GNU compiler.
        ! 101 is the number that fixes this issue for both compilers.
        

        block
            real(RK) :: unifrnd(101)
            call random_number(unifrnd)
!block
!integer(IK), allocatable :: RandomSeedValue(:)
!allocate(RandomSeedValue(RandomSeed%size))
!call random_seed(get=RandomSeedValue)
!write(*,"(*(g0,:,' '))") "unifrnd", unifrnd, RandomSeedValue
!end block
!if (this_image()==1) then
!    write(*,*) "RandomSeedObj%imageID, unifrnd: ", unifrnd
!    sync images(*)
!else
!    sync images(1)
!    write(*,*) "RandomSeedObj%imageID, unifrnd: ", unifrnd
!end if
!if (this_image()==1) read(*,*)
!sync all
        end block

        !else
            !call random_init( repeatable = RandomSeed%isRepeatable &
            !                , image_distinct = RandomSeed%isImageDistinct &
            !                , info = RandomSeed%info &
            !                , Err = RandomSeed%Err &
            !                , ProcessID = RandomSeed%ProcessID &
            !                )
        !end if

    end subroutine setRandomSeed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    ! This subroutine is not used anymore
!    ! This subroutine must be called by all images of a team
!    subroutine random_init(repeatable, image_distinct, info, Err, ProcessID)
!
!        use iso_fortran_env, only: int64
!#if defined IFORT_ENABLED
!        use ifport
!#endif
!        use Err_mod, only: Err_type
!        use Constants_mod, only: IK
!
!        implicit none
!
!        logical, intent(in), optional                       :: repeatable, image_distinct
!        character(:), allocatable, intent(out), optional    :: info
!        type(Err_type), intent(out), optional               :: Err
!        integer(IK) , intent(out), optional                 :: ProcessID
!
!        character(*), parameter                             :: PROCEDURE_NAME = "@random_init()"
!
!        logical                     :: isRepeatable, isImageDistinct, errIsPresent
!        integer(IK), allocatable    :: SeedValue(:)
!        integer(IK)                 :: i, seedSize, DateTimeValues(8) ! , iostat, fileUnit
!        integer(IK)                 :: pid = -huge(0)
!#if defined CAF_ENABLED
!        integer(IK)   , save        :: co_pid[*]
!        integer(int64), save        :: co_time[*] = -huge(0)
!#else
!        integer(IK)   , save        :: co_pid
!        integer(int64), save        :: co_time = -huge(0)
!#endif
!        integer(int64)              :: lcgInput
!
!        errIsPresent = present(Err)
!        if (errIsPresent) then
!            Err%occurred = .false.
!            Err%msg = ""
!        end if
!
!        isRepeatable = .true.
!        if (present(repeatable)) isRepeatable = repeatable
!
!        isImageDistinct = .false.
!        if (present(image_distinct)) isImageDistinct = image_distinct
!
!        call random_seed(size = seedSize)
!        allocate(SeedValue(seedSize))
!
!        if (isRepeatable) then
!            if (pid==-huge(0)) pid = getpid()
!        else
!            pid = getpid()
!        end if
!
!        if (present(ProcessID)) ProcessID = pid
!
!        ! First try if the OS provides a random number generator
!        !open( newunit   =   fileUnit &
!        !    , file      =   "/dev/urandom" &
!        !    , access    =   "stream" &
!        !    , form      =   "unformatted" &
!        !    , action    =   "read" &
!        !    , status    =   "old" &
!        !    , iostat    =   iostat &
!        !    )
!        !
!        !if (iostat == 0) then
!        !
!        !    if (present(info)) info = "OS provides random number generator."
!        !
!        !    if (errIsPresent) then
!        !        read(fileUnit,iostat=Err%stat) SeedValue
!        !        if (Err%stat/=0) then
!        !            Err%occurred = .true.
!        !            Err%msg = PROCEDURE_NAME // "Error occurred while reading array SeedValue from file='/dev/urandom'."
!        !            return
!        !        end if
!        !        close(fileUnit,iostat=Err%stat)
!        !        if (Err%stat/=0) then
!        !            Err%occurred = .true.
!        !            Err%msg = PROCEDURE_NAME // "Error occurred while attempting to close file='/dev/urandom'."
!        !            return
!        !        end if
!        !    else
!        !        read(fileUnit) SeedValue
!        !        close(fileUnit)
!        !    end if
!        !
!        !else
!
!            if (present(info)) info = "Ignoring the OS random number generator."
!
!            ! Fallback to XOR:ing the current time and co_pid. The co_pid is
!            ! useful in case one launches multiple instances of the same program in parallel.
!            if ( isImageDistinct ) then
!                if (isRepeatable) then
!                    if (co_time==-huge(0)) call getTime()
!                else
!                    call getTime()
!                end if
!            else
!#if defined CAF_ENABLED
!                if (this_image()==1) then
!#endif
!                    if (isRepeatable) then
!                        if (co_time==-huge(0)) call getTime()
!                    else
!                        call getTime()
!                    end if
!#if defined CAF_ENABLED
!                    sync images(*)
!                else
!                    sync images(1)
!                    co_time = co_time[1]
!                end if
!#endif
!            end if
!
!            if ( isImageDistinct ) then
!                co_pid = pid
!            else
!#if defined CAF_ENABLED
!                if (this_image()==1) then
!                    co_pid = pid
!                    sync images(*)
!                else
!                    sync images(1)
!                    co_pid = co_pid[1]
!                end if
!#else
!                co_pid = pid
!#endif
!            end if
!
!            lcgInput = ieor(co_time, int(co_pid, kind(co_time)))
!            do i = 1, seedSize
!                SeedValue(i) = lcg(lcgInput)
!            end do
!
!        !end if
!
!        call random_seed(put=SeedValue)
!
!    contains
!
!        ! This simple PRNG might not be good enough for real work, but is
!        ! sufficient for seeding a better PRNG.
!        function lcg(s)
!            integer :: lcg
!            integer(int64) :: s
!            if (s == 0) then
!            s = 104729
!            else
!            s = mod(s, 4294967296_int64)
!            end if
!            s = mod(s * 279470273_int64, 4294967291_int64)
!            lcg = int(mod(s, int(huge(0), int64)), kind(0))
!        end function lcg
!
!        subroutine getTime()
!            implicit none
!            call system_clock( count=co_time )
!            if (co_time <= 0) then
!                call date_and_time(values=DateTimeValues)
!                co_time = (DateTimeValues(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
!                        +  DateTimeValues(2) * 31_int64 * 24 * 60 * 60 * 1000 &
!                        +  DateTimeValues(3) * 24_int64 * 60 * 60 * 1000 &
!                        +  DateTimeValues(5) * 60 * 60 * 1000 &
!                        +  DateTimeValues(6) * 60 * 1000 &
!                        +  DateTimeValues(7) * 1000 &
!                        +  DateTimeValues(8)
!            end if
!        end subroutine getTime
!
!    end subroutine random_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module RandomSeed_mod