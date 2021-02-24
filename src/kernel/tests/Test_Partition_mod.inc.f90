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

!>  \brief This file implements the bodies of the partition modules.
!>  \author Amir Shahmoradi

    use Test_mod, only: Test_type

#if defined MAXDEN
    use PartitionMaxDen_mod
#elif defined MINVOL
    use PartitionMinVol_mod
#endif

    implicit none

    private
    public :: test_Partition

    type(Test_type) :: Test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 1000
        integer(IK)             :: nemax
        real(RK), allocatable   :: Domain(:,:)
        real(RK), allocatable   :: DomainSize(:)
        real(RK), allocatable   :: Point(:,:)
    contains
        procedure, pass :: read => readTestData
    end type TestData_type

    type(TestData_type) :: TestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Partition()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call TestData%read()
        call Test%run(test_runPartition_1, "test_runPartition_1")
        call Test%run(test_runPartition_2, "test_runPartition_2")
        !call Test%run(test_runPartition_3, "test_runPartition_3")
        !call Test%run(test_runPartition_4, "test_runPartition_4")
        !call Test%run(test_benchmark_1, "test_benchmark_1")
        !call Test%run(test_optimizePartition_1, "test_optimizePartition_1")
        call Test%finalize()
    end subroutine test_Partition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Generate a uniformly random set of points from an nd-cube.
    subroutine readTestData(TestData)
        implicit none
        class(TestData_type), intent(inout) :: TestData
        integer(IK) :: ip
        if (allocated(TestData%Point)) deallocate(TestData%Point) ! LCOV_EXCL_LINE
        if (allocated(TestData%Domain)) deallocate(TestData%Domain) ! LCOV_EXCL_LINE
        allocate(TestData%Point(TestData%nd,TestData%np))
        allocate(TestData%Domain(TestData%nd,2))
        TestData%Domain(:,1) = [0._RK, 0._RK]
        TestData%Domain(:,2) = [1._RK, 2._RK]
        TestData%DomainSize = TestData%Domain(:,2) - TestData%Domain(:,1)
        call random_number(TestData%Point)
        TestData%nemax = TestData%np / (TestData%nd + 1)
        do ip = 1, TestData%np
            TestData%Point(:,ip) = TestData%Domain(:,1) + TestData%Point(:,ip) * TestData%DomainSize
        end do
    end subroutine readTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runPartition()` by passing a fixed initial set of cluster centers to the Partition constructor.
    function test_runPartition_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        type(Partition_type)        :: Partition
        real(RK), allocatable       :: NormedPoint(:)
        real(RK)                    :: mahalSq
        integer                     :: ip, ic
        logical                     :: isInside

        assertion = .true.

        Partition = Partition_type  ( Point = TestData%Point & ! LCOV_EXCL_LINE
                                    , nsim = 0_IK & ! LCOV_EXCL_LINE
                                   !, nemax = TestData%nemax & ! LCOV_EXCL_LINE
                                    , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                    , logExpansion = log(1._RK) & ! LCOV_EXCL_LINE
                                   !, inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                    , parentLogVolNormed = sum(log(TestData%DomainSize)) & ! LCOV_EXCL_LINE
                                    , MaxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                    , MaxAllowedKmeansRecursion = 10000 & ! LCOV_EXCL_LINE
                                   !, maxAllowedKvolumeRecursion = 1000 & ! LCOV_EXCL_LINE
                                    )

        ! write data to output for further investigation

        Test%File = Test%openFile(label = "Partition")
        call Partition%write(Test%File%unit, TestData%nd, TestData%np, TestData%Point)
        close(Test%File%unit)

        assertion = assertion .and. .not. Partition%Err%occurred
        assertion = assertion .and. Partition%Err%stat /= 1_IK
        assertion = assertion .and. Partition%Err%stat /= 2_IK
        assertion = assertion .and. all(Partition%Size > 0_IK)
        assertion = assertion .and. all(Partition%Membership > 0_IK) .and. all(Partition%Membership <= Partition%neopt)

        do ip = 1, Partition%np
            ic = Partition%Membership(ip)
            if (Test%isVerboseMode .and. ic == 0_IK) write(Test%outputUnit,"(*(g0.15,:,', '))") "ic, nc, ip, np = ", ic, Partition%nc, ip, Partition%np
            NormedPoint = TestData%Point(:,ip) - Partition%Center(:,ic)
            mahalSq = dot_product(NormedPoint,matmul(Partition%InvCovMat(:,:,ic),NormedPoint))
            isInside = mahalSq - 1._RK <= 1.e-6_RK
            if (.not. isInside) then
                if (Test%isVerboseMode) write(Test%outputUnit,"(*(g0.15,:,' '))") "Point not inside!, mahalSq = ", mahalSq
                assertion = assertion .and. isInside
                exit
            end if
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%neopt                =", Partition%neopt
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Size < 1             =", pack(Partition%Size, mask = Partition%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Membership < 1       =", pack(Partition%Membership, mask = Partition%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Membership > neopt   =", pack(Partition%Membership, mask = Partition%Membership > Partition%neopt)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%occurred         =", Partition%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%stat             =", Partition%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%msg              =", Partition%Err%msg
            if (.not. isInside) write(Test%outputUnit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runPartition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runPartition()` by passing a fixed initial set of cluster centers to the Partition constructor.
    function test_runPartition_2() result(assertion)
        use iso_fortran_env, only: compiler_version
        use ClusteredPoint_mod, only: ClusteredPoint_type
        use RandomSeed_mod, only: RandomSeed_type
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        type(Partition_type)        :: Partition
        type(ClusteredPoint_type)   :: ClusteredPoint
        type(RandomSeed_type)       :: RandomSeed
        real(RK), allocatable       :: NormedPoint(:)
        real(RK)                    :: mahalSq
        logical                     :: isInside
        integer                     :: ip, ic

        integer(IK)                 :: itest, ntest
        namelist /specTest/ ntest

        character(:), allocatable   :: dist
        logical                     :: isRepeatable
        integer(IK)                 :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                    :: etamin, etamax, centerMin, centerMax
        namelist /specData/ rngseed, isRepeatable, nd, nc, sizeMin, sizeMax, etamin, etamax, centerMin, centerMax, dist

        integer(IK)                 :: nt, nsim, nemax, minSize
        integer(IK)                 :: maxAllowedKvolumeRecursion
        real(RK)                    :: mahalSqWeightExponent
        real(RK)                    :: inclusionFraction
        real(RK)                    :: expansionMaxDen
        real(RK)                    :: expansionMinVol
        logical                     :: rinitEnabled
        logical                     :: stanEnabled
        namelist /specPartition/ rngseed, isRepeatable, nc, nt, nemax, nsim, minSize, inclusionFraction, stanEnabled
        namelist /specPartition/ maxAllowedKvolumeRecursion, expansionMaxDen, expansionMinVol, mahalSqWeightExponent

        assertion = .true.

        rinitEnabled = compiler_version() == "GCC version 10.2.0"

        ! read the number of tests

        ntest = 1
        open(newunit = Test%File%unit, file = Test%inDir//"/Test_Partition_mod@test_runPartition_2.nml", status = "old")
        read(Test%File%unit, nml = specTest)
        close(Test%File%unit)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do itest = 1, ntest

            ! read cluster spec

            isRepeatable = .false.
            rngseed = -huge(rngseed)
            if (.not. allocated(dist)) allocate(character(63) :: dist)
            open(newunit = Test%File%unit, file = Test%inDir//"/Test_Partition_mod@test_runPartition_2.nml", status = "old")
            read(Test%File%unit, nml = specData)
            dist = trim(adjustl(dist))
            close(Test%File%unit)
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true., image_distinct = .true.)
            else
                if(rngseed /= -huge(rngseed)) RandomSeed = RandomSeed_type(imageID = Test%Image%id, inputSeed = rngseed)
            end if

            call ClusteredPoint%get ( nd = nd & ! LCOV_EXCL_LINE
                                    , nc = nc & ! LCOV_EXCL_LINE
                                    , etamin = etamin & ! LCOV_EXCL_LINE
                                    , etamax = etamax & ! LCOV_EXCL_LINE
                                    , sizeMin = sizeMin & ! LCOV_EXCL_LINE
                                    , sizeMax = sizeMax & ! LCOV_EXCL_LINE
                                    , centerMin = centerMin & ! LCOV_EXCL_LINE
                                    , centerMax = centerMax & ! LCOV_EXCL_LINE
                                    !, Size = [50, 1000, 500, 2000, 3] & ! LCOV_EXCL_LINE
                                    !, Eta = [1._RK, 2._RK, 0.5_RK, 0.05_RK, 1.5_RK] & ! LCOV_EXCL_LINE
                                    , dist = dist & ! LCOV_EXCL_LINE
                                    )

            ! write data to output for further investigation

            Test%File = Test%openFile(label = "ClusteredPoint")
            call ClusteredPoint%write(Test%File%unit)
            close(Test%File%unit)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Partition the clustered points clustered points
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! read partition spec

            nsim = 0_IK
            minSize = nd + 1_IK
            stanEnabled = .true.
            isRepeatable = .false.
            rngseed = -huge(rngseed)
            mahalSqWeightExponent = 0._RK
            maxAllowedKvolumeRecursion = 100_IK
            nemax = ClusteredPoint%np / (ClusteredPoint%nd + 1)
            open(newunit = Test%File%unit, file = Test%inDir//"/Test_Partition_mod@test_runPartition_2.nml", status = "old")
            read(Test%File%unit, nml = specPartition)
            close(Test%File%unit)
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true., image_distinct = .true.)
            else
                if(rngseed /= -huge(rngseed)) RandomSeed = RandomSeed_type(imageID = Test%Image%id, inputSeed = rngseed)
            end if

            !write(*,*) "zeroth nemax", nemax, ClusteredPoint%np, ClusteredPoint%nd + 1
            Partition = Partition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                        , nc = nc & ! LCOV_EXCL_LINE
                                        , nt = nt & ! LCOV_EXCL_LINE
                                        , nsim = nsim & ! LCOV_EXCL_LINE
                                        , nemax = nemax & ! LCOV_EXCL_LINE
                                        , trimEnabled = .false. & ! LCOV_EXCL_LINE
!#if defined MINVOL
!                                        , stanEnabled = .false. & ! LCOV_EXCL_LINE
!#endif
                                        , stanEnabled = stanEnabled & ! LCOV_EXCL_LINE
#if defined MINVOL
                                        , logExpansion = log(expansionMinVol) & ! LCOV_EXCL_LINE
#elif defined MAXDEN
                                        , logExpansion = log(expansionMaxDen) & ! LCOV_EXCL_LINE
#endif
                                        , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                        , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! LCOV_EXCL_LINE
                                       !, nemax = ClusteredPoint%nemax & ! LCOV_EXCL_LINE
                                       !, maxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                       !, maxAllowedKmeansRecursion = 10000 & ! LCOV_EXCL_LINE
                                       !, mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
                                        , maxAllowedKvolumeRecursion = maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                        , minSize = minSize & ! LCOV_EXCL_LINE
                                        )

            ! write data to output for further investigation

            Test%File = Test%openFile(label = "Partition")
            call Partition%write(Test%File%unit, ClusteredPoint%nd, ClusteredPoint%np, ClusteredPoint%Point)
            close(Test%File%unit)

            assertion = assertion .and. .not. Partition%Err%occurred
            assertion = assertion .and. Partition%Err%stat /= 1_IK
            assertion = assertion .and. Partition%Err%stat /= 2_IK
            assertion = assertion .and. all(Partition%Size(1:Partition%neopt) >= 0_IK)
            assertion = assertion .and. all(Partition%Membership > 0_IK) .and. all(Partition%Membership < Partition%neopt + 1)

            !if (.not. assertion) return

            ! Ensure all points are within their corresponding clusters.

            do ip = 1, Partition%np
                ic = Partition%Membership(ip)
                NormedPoint = ClusteredPoint%Point(:,ip) - Partition%Center(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(Partition%InvCovMat(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    if (Test%isVerboseMode) write(Test%outputUnit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                    assertion = assertion .and. isInside
                    exit
                end if
            end do

            if (Test%isVerboseMode .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%neopt                =", Partition%neopt
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Size < 1             =", pack(Partition%Size(1:Partition%neopt), mask = Partition%Size(1:Partition%neopt) < 1_IK)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Membership < 1       =", pack(Partition%Membership, mask = Partition%Membership < 1_IK)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Membership > neopt   =", pack(Partition%Membership, mask = Partition%Membership > Partition%neopt)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%occurred         =", Partition%Err%occurred
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%stat             =", Partition%Err%stat
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%Err%msg              =", Partition%Err%msg
                if (.not. isInside) write(Test%outputUnit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if

            if (.not. assertion) return

            ! The membership IDs must cover a full range from 1 to Partition%neopt

            block

                use Unique_mod, only: findUnique
                integer(IK) :: lenUnique
                integer(IK), allocatable :: UniqueValue(:), UniqueCount(:)

                call findUnique ( lenVector = Partition%np & ! LCOV_EXCL_LINE
                                , Vector = Partition%Membership & ! LCOV_EXCL_LINE
                                , lenUnique = lenUnique & ! LCOV_EXCL_LINE
                                , UniqueValue = UniqueValue & ! LCOV_EXCL_LINE
                                , UniqueCount = UniqueCount & ! LCOV_EXCL_LINE
                                )

                do ic = 1, Partition%neopt
                    assertion = assertion .and. any(UniqueValue == ic)
                    if (Test%isVerboseMode .and. .not. assertion) then
                        ! LCOV_EXCL_START
                        write(Test%outputUnit,"(*(g0.15,:,' '))")
                        write(Test%outputUnit,"(*(g0.15,:,' '))") "Membership IDs must cover a full range 1:Partition%neopt"
                        write(Test%outputUnit,"(*(g0.15,:,' '))") "Partition%neopt          =", Partition%neopt
                        write(Test%outputUnit,"(*(g0.15,:,' '))") "Unique Membership IDs    =", UniqueValue
                        write(Test%outputUnit,"(*(g0.15,:,' '))")
                        ! LCOV_EXCL_STOP
                    end if
                    if (.not. assertion) return
                end do

            end block

            if (.not. assertion) return

        end do

    end function test_runPartition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
