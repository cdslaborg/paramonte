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

module Test_RandomSeed_mod

    use RandomSeed_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_RandomSeed

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_RandomSeed()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call test_RandomSeed_type()
        call Test%finalize()

end subroutine test_RandomSeed

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_RandomSeed_type()

        implicit none
        type(RandomSeed_type)       :: RandomSeed
        integer                     :: seedSize
#if defined CAF_ENABLED
        integer, allocatable, save  :: Seed(:)[:]
#else
        integer, allocatable, save  :: Seed(:)
#endif



        if (Test%Image%isFirst) call Test%testing(   "RandomSeed_type with no input arguments to constructor &
                                                    &(default non-repeatable simulation, image-distinct)")

        RandomSeed = RandomSeed_type(imageID=Test%Image%id)  ! without any arguments it must give non-repeatable image-distinct random seeds.
        call Test%checkForErr(RandomSeed%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isRepeatable     : ", RandomSeed%isRepeatable
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isImageDistinct  : ", RandomSeed%isImageDistinct
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%info             : ", RandomSeed%info
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%size             : ", RandomSeed%size
            write(Test%outputUnit,"(*(g0))")
        end if

        call random_seed(size=seedSize)
        Test%assertion = RandomSeed%size == seedSize
        call Test%verify()

        if (Test%isDebugMode) then
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#if defined CAF_ENABLED
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#endif
            end if
#if defined CAF_ENABLED
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if

#if defined CAF_ENABLED
        allocate( Seed(seedSize)[*] )
#else
        allocate( Seed(seedSize) )
#endif

        call random_seed(get=Seed)
        Test%assertion = all(RandomSeed%Value == Seed)
        call Test%verify()




        if (Test%Image%isFirst) call Test%testing("RandomSeed_type with default isRepeatable, isImageDistinct=.false.")

        RandomSeed = RandomSeed_type(imageID=Test%Image%id, isImageDistinct=.false.)
        call Test%checkForErr(RandomSeed%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isRepeatable     : ", RandomSeed%isRepeatable
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isImageDistinct  : ", RandomSeed%isImageDistinct
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%info             : ", RandomSeed%info
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%size             : ", RandomSeed%size
            write(Test%outputUnit,"(*(g0))")
        end if

        if (Test%isDebugMode) then
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#if defined CAF_ENABLED
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#endif
            end if
#if defined CAF_ENABLED
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if

        call random_seed(get=Seed)
        Test%assertion = all(RandomSeed%Value == Seed)
        call Test%verify()


        if (Test%Image%isFirst) call Test%testing("RandomSeed_type for equivalence of Seed vector on all images")
        Test%assertion = .true.
#if defined CAF_ENABLED
        sync all
        if (Test%Image%id==1) then
            sync images(*)
        else
            if ( any(Seed /= Seed(:)[1]) ) Test%assertion = .true.
            sync images(1)
        end if
#endif
        call Test%verify()

#if defined CAF_ENABLED
        sync all
#endif




        if (Test%Image%isFirst) call Test%testing("RandomSeed_type with isRepeatable=.true., isImageDistinct=.false.")

        RandomSeed = RandomSeed_type(imageID=Test%Image%id, isRepeatable=.true., isImageDistinct=.false.)
        call Test%checkForErr(RandomSeed%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isRepeatable     : ", RandomSeed%isRepeatable
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isImageDistinct  : ", RandomSeed%isImageDistinct
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%info             : ", RandomSeed%info
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%size             : ", RandomSeed%size
            write(Test%outputUnit,"(*(g0))")
        end if

        if (Test%isDebugMode) then
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#if defined CAF_ENABLED
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#endif
            end if
#if defined CAF_ENABLED
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if

        call random_seed(get=Seed)
        Test%assertion = all(RandomSeed%Value == Seed)
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("RandomSeed_type for equivalence of Seed vector on all images")
        Test%assertion = .true.
#if defined CAF_ENABLED
        sync all
        if (Test%Image%id==1) then
            sync images(*)
        else
            if ( any(Seed /= Seed(:)[1]) ) Test%assertion = .false.
            sync images(1)
        end if
#endif
        call Test%verify()

        RandomSeed = RandomSeed_type(imageID=Test%Image%id, isRepeatable=.true., isImageDistinct=.false.)
        call Test%checkForErr(RandomSeed%Err)
        block
            integer, allocatable, save  :: SeedNew(:)
            allocate( SeedNew(seedSize) )
            call random_seed(get=SeedNew)
            if (Test%Image%isFirst) call Test%testing(   "RandomSeed_type for equivalence of the old and the new Seed vector on each image")
            if (Test%isDebugMode) then
                if (Test%Image%id==1) then
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedOld: ", Test%Image%id,",", Seed
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedNew: ", Test%Image%id,",", SeedNew
#if defined CAF_ENABLED
                else
                    sync images (Test%Image%id-1)
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedOld: ", Test%Image%id,",", Seed
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedNew: ", Test%Image%id,",", SeedNew
#endif
                end if
#if defined CAF_ENABLED
                if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
            end if
            Test%assertion = all(SeedNew==Seed)
            call Test%verify()
            deallocate(SeedNew)
        end block

#if defined CAF_ENABLED
        sync all
#endif



        if (Test%Image%isFirst) call Test%testing("RandomSeed_type with isRepeatable=.true., isImageDistinct=.true.")

        RandomSeed = RandomSeed_type(imageID=Test%Image%id, isRepeatable=.true., isImageDistinct=.true.)
        call Test%checkForErr(RandomSeed%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isRepeatable     : ", RandomSeed%isRepeatable
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%isImageDistinct  : ", RandomSeed%isImageDistinct
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%info             : ", RandomSeed%info
            write(Test%outputUnit,"(*(g0))")        "RandomSeed%size             : ", RandomSeed%size
            write(Test%outputUnit,"(*(g0))")
        end if

        if (Test%isDebugMode) then
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#if defined CAF_ENABLED
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), GETPID(), Seed: ", Test%Image%id, ",", RandomSeed%imageID, ",", RandomSeed%Value
#endif
            end if
#if defined CAF_ENABLED
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if

        call random_seed(get=Seed)
        Test%assertion = all(RandomSeed%Value == Seed)
        call Test%verify()
#if defined CAF_ENABLED
        sync all
#endif

        RandomSeed = RandomSeed_type(imageID=Test%Image%id, isRepeatable=.true., isImageDistinct=.true.)
        call Test%checkForErr(RandomSeed%Err)
        block
            integer, allocatable, save  :: SeedNew(:)
            allocate( SeedNew(seedSize) )
            call random_seed(get=SeedNew)
            if (Test%Image%isFirst) call Test%testing("RandomSeed_type for equivalence of the old and the new Seed vector on each image")
            if (Test%isDebugMode) then
                if (Test%Image%id==1) then
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedOld(diff. on each image): ", Test%Image%id,",", Seed
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedNew(diff. on each image): ", Test%Image%id,",", SeedNew
#if defined CAF_ENABLED
                else
                    sync images (Test%Image%id-1)
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedOld(diff. on each image): ", Test%Image%id,",", Seed
                    write(Test%outputUnit,"(*(g0,' '))") "this_image(), SeedNew(diff. on each image): ", Test%Image%id,",", SeedNew
#endif
                end if
#if defined CAF_ENABLED
                if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
            end if
            Test%assertion = all(SeedNew==Seed)
            call Test%verify()
            deallocate(SeedNew)
        end block

        if (Test%Image%isFirst) call Test%testing("RandomSeed_type for non-equivalence of Seed vector on all images")
        Test%assertion = .true.
#if defined CAF_ENABLED
        sync all
        if (Test%Image%id==1) then
            sync images(*)
        else
            if ( all(Seed == Seed(:)[1]) ) Test%assertion = .false.
            sync images(1)
        end if
#endif
        call Test%verify()




        if (Test%Image%isFirst) call Test%testing("RandomSeed_type(inputSeed = 1313, isImageDistinct=.false.)")
        RandomSeed = RandomSeed_type(imageID=Test%Image%id, inputSeed = 1313, isImageDistinct=.false.)
        call Test%checkForErr(RandomSeed%Err)
        if (Test%isDebugMode) then
#if defined CAF_ENABLED
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), Seed(same on each image): ", Test%Image%id,",", RandomSeed%Value
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), Seed(same on each image): ", Test%Image%id,",", RandomSeed%Value
            end if
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if
        call random_seed(get=Seed)
        Test%assertion = .true.
#if defined CAF_ENABLED
        sync all
        if (Test%Image%id==1) then
            sync images(*)
        else
            if ( any(Seed /= Seed(:)[1]) ) Test%assertion = .false.
            sync images(1)
        end if
#endif
        call Test%verify()




        if (Test%Image%isFirst) call Test%testing("RandomSeed_type(inputSeed = 1313, isImageDistinct=.true.)")
        RandomSeed = RandomSeed_type(imageID=Test%Image%id, inputSeed = 1313, isImageDistinct=.true.)
        call Test%checkForErr(RandomSeed%Err)
        if (Test%isDebugMode) then
            if (Test%Image%id==1) then
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), Seed(same on each image): ", Test%Image%id,",", RandomSeed%Value
#if defined CAF_ENABLED
            else
                sync images (Test%Image%id-1)
                write(Test%outputUnit,"(*(g0,' '))") "this_image(), Seed(same on each image): ", Test%Image%id,",", RandomSeed%Value
#endif
            end if
#if defined CAF_ENABLED
            if (Test%Image%id<Test%Image%count) sync images (Test%Image%id+1)
#endif
        end if
        call random_seed(get=Seed)
        Test%assertion = .true.
#if defined CAF_ENABLED
        sync all
        if (Test%Image%id==1) then
            sync images(*)
        else
            if ( any(Seed == Seed(:)[1]) ) Test%assertion = .false.
            sync images(1)
        end if
#endif
        call Test%verify()


        deallocate( Seed )

    end subroutine test_RandomSeed_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_RandomSeed_mod