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

module SpecBase_RandomSeed_mod

    use RandomSeed_mod, only: RandomSeed_t => RandomSeed_type
    use Constants_mod, only: IK
    implicit none

#if defined CAF_ENABLED
    type(RandomSeed_t), save        :: comv_RandomSeed[*]
#else
    type(RandomSeed_t), save        :: comv_RandomSeed
#endif

    character(*), parameter :: MODULE_NAME = "@SpecBase_RandomSeed_mod"

    integer(IK)                     :: randomSeed ! namelist input

    type                            :: RandomSeed_type
        logical                     :: isImageDistinct
        logical                     :: isRepeatable
        integer(IK)                 :: userSeed
        integer(IK)                 :: nullSeed
        integer(IK)                 :: sizeSeed
        integer(IK)                 :: imageID, imageCount
        integer(IK)                 :: ProcessID
        integer(IK) , allocatable   :: Seed(:,:)
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setRandomSeed, nullifyNameListVar   ! , checkForSanity
    end type RandomSeed_type

    interface RandomSeed_type
        module procedure            :: constructRandomSeed
    end interface RandomSeed_type
    private :: constructRandomSeed, setRandomSeed, nullifyNameListVar   ! , checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructRandomSeed(methodName,imageID,imageCount) result(RandomSeedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructRandomSeed
#endif
        use Constants_mod, only: IK, NULL_IK
        implicit none
        character(*), intent(in)    :: methodName
        integer(IK), intent(in)     :: imageID, imageCount
        type(RandomSeed_type)       :: RandomSeedObj
        RandomSeedObj%userSeed          = NULL_IK
        RandomSeedObj%nullSeed          = NULL_IK
        RandomSeedObj%ProcessID         = NULL_IK
        RandomSeedObj%isRepeatable      = .false.
        RandomSeedObj%isImageDistinct   = .true.
        RandomSeedObj%imageID           = imageID
        RandomSeedObj%imageCount        = imageCount
        call random_seed(size=RandomSeedObj%sizeSeed)
        allocate(RandomSeedObj%Seed(RandomSeedObj%sizeSeed,RandomSeedObj%imageCount))
        RandomSeedObj%desc = &
        "randomSeed is a scalar 32bit integer that serves as the seed of the random number generator. When it is provided, &
        &the seed of the random number generator will be set in a specific deterministic manner to enable future replications &
        &of the simulation with the same configuration and input specifications. The default value for randomSeed is an integer &
        &vector of processor-dependent size and value that will vary from one simulation to another. &
        &However, enough care has been taken to assign unique random seed values to the random number generator on &
        &each of the parallel threads (or images, processors, cores, ...) at all circumstances."
    end function constructRandomSeed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(RandomSeedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(RandomSeed_type), intent(in) :: RandomSeedObj
        randomSeed = RandomSeedObj%nullSeed
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setRandomSeed(RandomSeedObj,randomSeed,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setRandomSeed
#endif

        use Constants_mod, only: IK
        use Err_mod, only: Err_type

        implicit none

        class(RandomSeed_type), intent(inout)   :: RandomSeedObj
        integer(IK), intent(in)                 :: randomSeed
        type(Err_type), intent(out)             :: Err

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@setRandomSeed()"

        RandomSeedObj%userSeed = randomSeed

        ! broadcast all random seed values of all images to all images

        blockBcastSeeds: block
#if defined CAF_ENABLED
            integer(IK) :: imageID
#elif defined MPI_ENABLED
            use mpi
            integer :: ierrMPI
            integer(IK), allocatable :: Seed(:,:)
#endif
            if ( RandomSeedObj%userSeed == RandomSeedObj%nullSeed ) then
                comv_RandomSeed = RandomSeed_t  ( imageID            = RandomSeedObj%imageID         &
                                                , isRepeatable       = RandomSeedObj%isRepeatable    &
                                                , isImageDistinct    = RandomSeedObj%isImageDistinct &
                                                )
            else
                comv_RandomSeed = RandomSeed_t  ( imageID            = RandomSeedObj%imageID         &
                                                , isRepeatable       = RandomSeedObj%isRepeatable    &
                                                , isImageDistinct    = RandomSeedObj%isImageDistinct &
                                                , inputSeed          = RandomSeedObj%userSeed        &
                                                )
            end if

#if defined GNU_ENABLED && CAF_ENABLED
            ! opencoarrays crashes without this, by somehow setting comv_RandomSeed%Err%occurred = TRUE
            ! likely a result of memory corruption
            if (comv_RandomSeed%Err%occurred) write(*,*) ""
#endif

            if (comv_RandomSeed%Err%occurred) then
                Err%occurred = .true.
                Err%msg = Err%msg // PROCEDURE_NAME // comv_RandomSeed%Err%msg
                return
            end if

            call comv_RandomSeed%get()
            RandomSeedObj%Seed(:,RandomSeedObj%imageID) = comv_RandomSeed%Value(:)
#if defined CAF_ENABLED
            sync all    ! allow all images to set the seed first, then fetch the values
            do imageID = 1, RandomSeedObj%imageCount
                if (imageID/=RandomSeedObj%imageID) RandomSeedObj%Seed(:,imageID) = comv_RandomSeed[imageID]%Value(:)
            end do
#elif defined MPI_ENABLED
            allocate(Seed(RandomSeedObj%sizeSeed,RandomSeedObj%imageCount))
            call mpi_barrier(mpi_comm_world,ierrMPI) ! allow all images to set the seed first, then fetch the values
            call mpi_allgather  ( RandomSeedObj%Seed(:,RandomSeedObj%imageID)   &   ! send buffer
                                , RandomSeedObj%sizeSeed                        &   ! send count
                                , mpi_integer                                   &   ! send datatype
                                , Seed(:,:)                                     &   ! receive buffer
                                , RandomSeedObj%sizeSeed                        &   ! receive count
                                , mpi_integer                                   &   ! receive datatype
                                , mpi_comm_world                                &   ! comm
                                , ierrMPI                                       &   ! ierr
                                )
            RandomSeedObj%Seed(:,:) = Seed
            deallocate(Seed)
#endif
        end block blockBcastSeeds

    end subroutine setRandomSeed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_RandomSeed_mod