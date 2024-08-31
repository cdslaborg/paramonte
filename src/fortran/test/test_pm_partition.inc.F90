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

!>  \brief This file implements the bodies of the partition modules.
!>  \author Amir Shahmoradi

    use pm_test, only: test_type, LK

#if defined MAXDEN_RECURSIVE
    use pm_partitionRecursiveMaxDen
#elif defined MINVOL_RECURSIVE
    use pm_partitionRecursiveMinVol
#elif defined MAXDEN
    use pm_partitionMaxDen
#elif BENCHM_ENABLED
    use Benchm_pmod, only: getPartition, MODULE_NAME
    use pm_partitionMaxDen, only: partition_type, RK
#endif

    implicit none

    private
    public :: setTest
    type(test_type) :: test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 1000
        integer(IK)             :: nemax
        real(RK), allocatable   :: Domain(:,:)
        real(RK), allocatable   :: DomainSize(:)
        real(RK), allocatable   :: Point(:,:)
    contains
        procedure, pass :: get => getTestData
    end type TestData_type

    type(TestData_type) :: TestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call TestData%get()
        call test%run(test_constructPartition_1, SK_"test_constructPartition_1")
#if defined MAXDEN
        call test%run(test_getVarCorrectionScaleFactorSq_1, SK_"test_getVarCorrectionScaleFactorSq_1")
#endif
#if !BENCHM_ENABLED
        !call test%run(test_constructPartition_2, SK_"test_constructPartition_2")
        !call test%run(test_constructPartition_3, SK_"test_constructPartition_3")
#endif
        !call test%run(test_constructPartition_4, SK_"test_constructPartition_4")
        !call test%run(test_benchmark_1, SK_"test_benchmark_1")
        !call test%run(test_optimizePartition_1, SK_"test_optimizePartition_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Generate a uniformly random set of points from an nd-cube.
    subroutine getTestData(TestData)
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
    end subroutine getTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `constructPartition()` by passing a fixed initial set of cluster centers to the Partition constructor.
    function test_constructPartition_1() result(assertion)
        use iso_fortran_env, only: compiler_version
        use ClusteredPoint_pmod, only: ClusteredPoint_type
        use pm_randomSeed, only: randomSeed_type
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
        type(partition_type)            :: Partition
        type(ClusteredPoint_type)       :: ClusteredPoint
        type(randomSeed_type)           :: randomSeed
        real(RK), allocatable           :: NormedPoint(:)
        real(RK)                        :: mahalSq
        logical(LK)                     :: isInside
        integer(IK)                     :: ip, ic
    
        integer(IK)                     :: itest, ntest
        namelist /specTest/ ntest

        character(:, SK), allocatable   :: dist
        logical(LK)                     :: isRepeatable
        integer(IK)                     :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                        :: etamin, etamax, centerMin, centerMax, relTol
        namelist /specData/ rngseed, isRepeatable, nd, nc, sizeMin, sizeMax, etamin, etamax, centerMin, centerMax, dist, relTol

        integer(IK)                     :: nt, nsim, nemax, minSize
        integer(IK)                     :: kvolumeNumRecursionMax
        real(RK)                        :: mahalSqWeightExponent
        real(RK)                        :: inclusionFraction
        real(RK)                        :: expansionMaxDen
        real(RK)                        :: expansionMinVol
        real(RK)                        :: shrinkageMaxDen
        real(RK)                        :: shrinkageMinVol
        logical(LK)                     :: rinitEnabled
        logical(LK)                     :: stanEnabled
        logical(LK)                     :: biasCorrectionEnabled
        logical(LK)                     :: shapeAdaptationEnabled, scaleOptimizationEnabled, shapeOptimizationEnabled
        namelist /specPartition/ rngseed, isRepeatable, nc, nt, nemax, nsim, minSize, inclusionFraction, stanEnabled
        namelist /specPartition/ expansionMaxDen, expansionMinVol, shrinkageMaxDen, shrinkageMinVol
        namelist /specPartition/ kvolumeNumRecursionMax, mahalSqWeightExponent
        namelist /specPartition/ biasCorrectionEnabled, shapeAdaptationEnabled, scaleOptimizationEnabled, shapeOptimizationEnabled

        assertion = .true._LK

        rinitEnabled = compiler_version() == "GCC version 10.2.0"

        ! read the number of tests

        ntest = 1
        open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
        read(test%File%unit, nml = specTest)
        close(test%File%unit)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do itest = 1, ntest

            ! read cluster spec

            relTol = 0.01_RK
            isRepeatable = .false._LK
            rngseed = -huge(rngseed)
            if (.not. allocated(dist)) allocate(character(63) :: dist)
            open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
            read(test%File%unit, nml = specData)
            dist = trim(adjustl(dist))
            close(test%File%unit)
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true._LK, image_distinct = .true._LK)
            else
                if(rngseed /= -huge(rngseed)) randomSeed = randomSeed_type(imageID = test%image%id, inputSeed = rngseed)
            end if
!nc = 4
!sizeMin = 1000
!sizeMax = 1000
!centerMin = -5._RK
!centerMax = 5._RK
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
                                    , hubEnabled = .false._LK & ! LCOV_EXCL_LINE
                                    , relTol = relTol & ! LCOV_EXCL_LINE
                                    )

            ! write data to output for further investigation

            test%File = test%openFile(label = "ClusteredPoint")
            call ClusteredPoint%write(test%File%unit)
            close(test%File%unit)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Partition the clustered points clustered points
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! read partition spec

            nsim = 0
            minSize = 0 !nd + 1
            stanEnabled = .true._LK
            isRepeatable = .false._LK
            rngseed = -huge(rngseed)
            expansionMinVol = 1._RK
            expansionMaxDen = 1._RK
            shrinkageMinVol = 1._RK
            shrinkageMaxDen = 1._RK
            mahalSqWeightExponent = 0._RK
            biasCorrectionEnabled = .true._LK
            shapeAdaptationEnabled = .true._LK
            scaleOptimizationEnabled = .true._LK
            shapeOptimizationEnabled = .true._LK
            kvolumeNumRecursionMax = 100
            nemax = ClusteredPoint%np / (ClusteredPoint%nd + 1)
            open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
            read(test%File%unit, nml = specPartition)
            close(test%File%unit)
            !if (shrinkageMaxDen > 0._RK) then; shrinkageMaxDen = log(shrinkageMaxDen); else; shrinkageMaxDen = NEGBIG_RK; endif ! LCOV_EXCL_LINE
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true._LK, image_distinct = .true._LK)
            else
                if(rngseed /= -huge(rngseed)) randomSeed = randomSeed_type(imageID = test%image%id, inputSeed = rngseed)
            end if

#if BENCHM_ENABLED

            Partition = getPartition( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                    , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! & ! LCOV_EXCL_LINE
                                    , minSize = minSize & ! LCOV_EXCL_LINE
                                    , nemax = nemax & ! LCOV_EXCL_LINE
                                    )

#else

            !write(*,*) "zeroth nemax", nemax, ClusteredPoint%np, ClusteredPoint%nd + 1
            Partition = partition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                        , nc = nc & ! LCOV_EXCL_LINE
                                        , nt = nt & ! LCOV_EXCL_LINE
                                        , nsim = nsim & ! LCOV_EXCL_LINE
                                        , nemax = nemax & ! LCOV_EXCL_LINE
                                        , trimEnabled = .false._LK & ! LCOV_EXCL_LINE
!#if defined MINVOL
!                                        , stanEnabled = .false._LK & ! LCOV_EXCL_LINE
!#endif
                                        , stanEnabled = stanEnabled & ! LCOV_EXCL_LINE
#if defined MINVOL_RECURSIVE
                                        , logExpansion = log(expansionMinVol) & ! LCOV_EXCL_LINE
                                        , logShrinkage = log(shrinkageMinVol) & ! LCOV_EXCL_LINE
                                        , kvolumeNumRecursionMax = min(kvolumeNumRecursionMax, 20) & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = min(mahalSqWeightExponent,1._RK) & ! LCOV_EXCL_LINE
#elif defined MAXDEN_RECURSIVE || defined MAXDEN
                                        , logExpansion = log(expansionMaxDen) & ! LCOV_EXCL_LINE
                                        , logShrinkage = log(shrinkageMaxDen) & ! LCOV_EXCL_LINE
                                        , scaleOptimizationEnabled = scaleOptimizationEnabled & ! LCOV_EXCL_LINE
                                        , shapeOptimizationEnabled = shapeOptimizationEnabled & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#if defined MAXDEN_RECURSIVE
                                        , kvolumeNumRecursionMax = min(kvolumeNumRecursionMax, 20) & ! LCOV_EXCL_LINE
                                        !, kvolumeNumRecursionMax = kvolumeNumRecursionMax & ! LCOV_EXCL_LINE
#else
                                        , kvolumeNumRecursionMax = kvolumeNumRecursionMax & ! LCOV_EXCL_LINE
                                        , shapeAdaptationEnabled = shapeAdaptationEnabled & ! LCOV_EXCL_LINE
                                        , biasCorrectionEnabled = biasCorrectionEnabled & ! LCOV_EXCL_LINE
#endif
#endif
#if defined MAXDEN_RECURSIVE || defined MINVOL_RECURSIVE
                                        , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
#endif
                                        , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! & ! LCOV_EXCL_LINE
                                       !, nemax = ClusteredPoint%nemax & ! LCOV_EXCL_LINE
                                       !, kmeansNumFailMax = 10000 & ! LCOV_EXCL_LINE
                                       !, kmeansNumRecursionMax = 10000 & ! LCOV_EXCL_LINE
                                       !, mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
                                        , minSize = minSize & ! LCOV_EXCL_LINE
                                        )

            ! write data to output for further investigation

            test%File = test%openFile(label = "Partition")
            call partition%write(test%File%unit, ClusteredPoint%Point)
            close(test%File%unit)

            assertion = assertion .and. .not. partition%err%occurred
            call test%assert(assertion)
            assertion = assertion .and. partition%err%stat /= 1
            call test%assert(assertion)
            assertion = assertion .and. partition%err%stat /= 2
            call test%assert(assertion)
            assertion = assertion .and. all(partition%Size(1:partition%ne) >= 0)
            call test%assert(assertion)
            assertion = assertion .and. all(partition%membership > 0) .and. all(partition%membership < partition%ne + 1)

            !if (.not. assertion) return

            ! Ensure all points are within their corresponding clusters.

            do ip = 1, partition%np
                ic = partition%membership(ip)
                NormedPoint = ClusteredPoint%Point(:,ip) - partition%Center(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(partition%invCov(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    if (test%traceable) write(test%disp%unit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                    assertion = assertion .and. isInside
                    exit
                end if
            end do
            call test%assert(assertion)

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne                 =", partition%ne
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%Size < 1           =", pack(partition%Size(1:partition%ne), mask = partition%Size(1:partition%ne) < 1)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership < 1     =", pack(partition%membership, mask = partition%membership < 1)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership > ne    =", pack(partition%membership, mask = partition%membership > partition%ne)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%occurred       =", partition%err%occurred
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%stat           =", partition%err%stat
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%msg            =", partition%err%msg
                if (.not. isInside) write(test%disp%unit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
                write(test%disp%unit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
            !assertion = .true._LK
            call test%assert(assertion)

            ! The membership IDs must cover a full range from 1 to partition%ne

            block

                use pm_arrayUnique, only: setUnique
                integer(IK) :: lenUnique
                integer(IK), allocatable :: UniqueValue(:), UniqueCount(:)

                call setUnique  ( Vector = partition%membership & ! LCOV_EXCL_LINE
                                , UniqueValue = UniqueValue & ! LCOV_EXCL_LINE
                                , UniqueCount = UniqueCount & ! LCOV_EXCL_LINE
                                )
                lenUnique = size(UniqueValue, kind = IK)
                do ic = 1, partition%ne
                    assertion = assertion .and. any(UniqueValue == ic)
                    if (test%traceable .and. .not. assertion) then
                        ! LCOV_EXCL_START
                        write(test%disp%unit,"(*(g0.15,:,' '))")
                        write(test%disp%unit,"(*(g0.15,:,' '))") "Membership IDs must cover a full range 1:partition%ne"
                        write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne             =", partition%ne
                        write(test%disp%unit,"(*(g0.15,:,' '))") "Unique Membership IDs    =", UniqueValue
                        write(test%disp%unit,"(*(g0.15,:,' '))")
                        ! LCOV_EXCL_STOP
                    end if
                    if (.not. assertion) exit
                end do
                call test%assert(assertion)

            end block

            if (.not. assertion) return
#endif

        end do

    end function test_constructPartition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `constructPartition()` by passing a fixed initial set of cluster centers to the Partition constructor.
    function test_constructPartition_2() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        use pm_ellipsoid, only: getLogVolUnitBall
        implicit none
        logical(LK)                 :: assertion
        type(partition_type)        :: Partition
        real(RK), allocatable       :: NormedPoint(:)
        real(RK)                    :: mahalSq
        integer(IK)                 :: ip, ic
        logical(LK)                 :: isInside

        assertion = .true._LK

        Partition = partition_type  ( Point = TestData%Point & ! LCOV_EXCL_LINE
                                    , nsim = 0 & ! LCOV_EXCL_LINE
                                    , nemax = TestData%np & ! LCOV_EXCL_LINE
                                    , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                   !, logExpansion = log(1._RK) & ! LCOV_EXCL_LINE
                                   !, inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                    , parentLogVolNormed = sum(log(TestData%DomainSize)) - getLogVolUnitBall(TestData%nd) & ! LCOV_EXCL_LINE
                                    , KmeansNumFailMax = 10000 & ! LCOV_EXCL_LINE
                                    , kmeansNumRecursionMax = 10000 & ! LCOV_EXCL_LINE
                                   !, kvolumeNumRecursionMax = 1000 & ! LCOV_EXCL_LINE
                                    )

        ! write data to output for further investigation

        test%File = test%openFile(label = "Partition")
        call partition%write(test%File%unit, TestData%Point)
        close(test%File%unit)

        assertion = assertion .and. .not. partition%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. partition%err%stat /= 1
        call test%assert(assertion)
        assertion = assertion .and. partition%err%stat /= 2
        call test%assert(assertion)
        assertion = assertion .and. all(partition%Size > 0)
        call test%assert(assertion)
        assertion = assertion .and. all(partition%membership > 0) .and. all(partition%membership <= partition%ne)
        call test%assert(assertion)

        do ip = 1, partition%np
            ic = partition%membership(ip)
            if (test%traceable .and. ic == 0) write(test%disp%unit,"(*(g0.15,:,', '))") "ic, nc, ip, np = ", ic, partition%nc, ip, partition%np
            NormedPoint = TestData%Point(:,ip) - partition%Center(:,ic)
            mahalSq = dot_product(NormedPoint,matmul(partition%invCov(:,:,ic),NormedPoint))
            isInside = mahalSq - 1._RK <= 1.e-6_RK
            if (.not. isInside) then
                if (test%traceable) write(test%disp%unit,"(*(g0.15,:,' '))") "Point not inside!, mahalSq = ", mahalSq
                assertion = assertion .and. isInside
                exit
            end if
        end do
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne                 =", partition%ne
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%Size < 1           =", pack(partition%Size, mask = partition%Size < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership < 1     =", pack(partition%membership, mask = partition%membership < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership > ne    =", pack(partition%membership, mask = partition%membership > partition%ne)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%occurred       =", partition%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%stat           =", partition%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%msg            =", partition%err%msg
            if (.not. isInside) write(test%disp%unit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructPartition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> For a very limited input `nemax`, `constructPartition()` must return `ne = nemax`.
    function test_constructPartition_3() result(assertion)
        use iso_fortran_env, only: compiler_version
        use ClusteredPoint_pmod, only: ClusteredPoint_type
        use pm_randomSeed, only: randomSeed_type
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        type(partition_type)        :: Partition
        type(ClusteredPoint_type)   :: ClusteredPoint
        type(randomSeed_type)       :: randomSeed
        real(RK), allocatable       :: NormedPoint(:)
        real(RK)                    :: mahalSq
        logical(LK)                 :: isInside
        integer(IK)                 :: ip, ic

        character(:, SK), allocatable   :: dist
        logical(LK)                 :: isRepeatable
        integer(IK)                 :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                    :: etamin, etamax, centerMin, centerMax

        integer(IK)                 :: nt, nsim, nemax, minSize
        integer(IK)                 :: kvolumeNumRecursionMax
        real(RK)                    :: mahalSqWeightExponent
        real(RK)                    :: inclusionFraction
        real(RK)                    :: expansionMaxDen
        real(RK)                    :: expansionMinVol
        real(RK)                    :: shrinkageMaxDen
        real(RK)                    :: shrinkageMinVol

        logical(LK)                 :: scaleOptimizationEnabled
        logical(LK)                 :: shapeOptimizationEnabled
        logical(LK)                 :: rinitEnabled
        logical(LK)                 :: stanEnabled

        assertion = .true._LK

        ! read the number of tests

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nd = 2
        nc = 10
        etamin = 0.5d0
        etamax = 10.d0
        sizeMin = 1000
        sizeMax = 100
        centerMin = -2.5d0
        centerMax = +2.5d0

        ! read cluster spec

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


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Partition the clustered points clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! read partition spec

        nc = getUnifRand(2,6)
        nemax = getUnifRand(1,10)

        Partition = partition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                    , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! LCOV_EXCL_LINE
                                    , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                    , nemax = nemax & ! LCOV_EXCL_LINE
                                    , nc = nc & ! LCOV_EXCL_LINE
                                    )

        assertion = assertion .and. .not. partition%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. partition%nemax >= partition%ne
        call test%assert(assertion)
        assertion = assertion .and. all(partition%Size(1:partition%ne) >= 0)
        call test%assert(assertion)
        assertion = assertion .and. all(partition%membership > 0) .and. all(partition%membership < partition%ne + 1)
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne                 =", partition%ne
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%nemax              =", partition%nemax
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%Size < 1           =", pack(partition%Size(1:partition%ne), mask = partition%Size(1:partition%ne) < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership < 1     =", pack(partition%membership, mask = partition%membership < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership > ne    =", pack(partition%membership, mask = partition%membership > partition%ne)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%occurred       =", partition%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%stat           =", partition%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%msg            =", partition%err%msg
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if

        ! Ensure all points are within their corresponding clusters.

        do ip = 1, partition%np
            ic = partition%membership(ip)
            NormedPoint = ClusteredPoint%Point(:,ip) - partition%Center(:,ic)
            mahalSq = dot_product(NormedPoint,matmul(partition%invCov(:,:,ic),NormedPoint))
            isInside = mahalSq - 1._RK <= 1.e-6_RK
            if (.not. isInside) then
                if (test%traceable) write(test%disp%unit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                assertion = assertion .and. isInside
                exit
            end if
        end do
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne                 =", partition%ne
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%Size < 1           =", pack(partition%Size(1:partition%ne), mask = partition%Size(1:partition%ne) < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership < 1     =", pack(partition%membership, mask = partition%membership < 1)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership > ne    =", pack(partition%membership, mask = partition%membership > partition%ne)
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%occurred       =", partition%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%stat           =", partition%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%msg            =", partition%err%msg
            if (.not. isInside) write(test%disp%unit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if

        ! The membership IDs must cover a full range from 1 to partition%ne

        block

            use pm_arrayUnique, only: setUnique
            integer(IK)                 :: lenUnique
            integer(IK) , allocatable   :: UniqueValue(:), UniqueCount(:)

            call setUnique  ( Array = partition%membership & ! LCOV_EXCL_LINE
                            , Unique = UniqueValue & ! LCOV_EXCL_LINE
                            , Count = UniqueCount & ! LCOV_EXCL_LINE
                            )
            lenUnique = size(UniqueValue, kind = IK)
            do ic = 1, partition%ne
                assertion = assertion .and. any(UniqueValue == ic)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0.15,:,' '))")
                    write(test%disp%unit,"(*(g0.15,:,' '))") "Membership IDs must cover a full range 1:partition%ne"
                    write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne             =", partition%ne
                    write(test%disp%unit,"(*(g0.15,:,' '))") "Unique Membership IDs    =", UniqueValue
                    write(test%disp%unit,"(*(g0.15,:,' '))")
                    ! LCOV_EXCL_STOP
                end if
                if (.not. assertion) return
            end do

        end block

        if (.not. assertion) return

    end function test_constructPartition_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test the performance of the methods for a uniform square in high dimensions.
    !> \todo this test is not yet implemented.
    function test_constructPartition_4() result(assertion)
        use iso_fortran_env, only: compiler_version
        use ClusteredPoint_pmod, only: ClusteredPoint_type
        use pm_randomSeed, only: randomSeed_type
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        type(partition_type)        :: Partition
        type(ClusteredPoint_type)   :: ClusteredPoint
        type(randomSeed_type)       :: randomSeed
        real(RK), allocatable       :: NormedPoint(:)
        real(RK)                    :: mahalSq
        logical(LK)                 :: isInside
        integer(IK)                 :: ip, ic

        integer(IK)                 :: itest, ntest
        namelist /specTest/ ntest

        character(:, SK), allocatable   :: dist
        logical(LK)                 :: isRepeatable
        integer(IK)                 :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                    :: etamin, etamax, centerMin, centerMax, relTol
        namelist /specData/ rngseed, isRepeatable, nd, nc, sizeMin, sizeMax, etamin, etamax, centerMin, centerMax, dist, relTol

        integer(IK)                 :: nt, nsim, nemax, minSize
        integer(IK)                 :: kvolumeNumRecursionMax
        real(RK)                    :: mahalSqWeightExponent
        real(RK)                    :: inclusionFraction
        real(RK)                    :: expansionMaxDen
        real(RK)                    :: expansionMinVol
        real(RK)                    :: shrinkageMaxDen
        real(RK)                    :: shrinkageMinVol
        logical(LK)                 :: rinitEnabled
        logical(LK)                 :: stanEnabled
        logical(LK)                 :: shapeAdaptationEnabled, scaleOptimizationEnabled, shapeOptimizationEnabled
        namelist /specPartition/ rngseed, isRepeatable, nc, nt, nemax, nsim, minSize, inclusionFraction, stanEnabled
        namelist /specPartition/ expansionMaxDen, expansionMinVol, shrinkageMaxDen, shrinkageMinVol
        namelist /specPartition/ kvolumeNumRecursionMax, mahalSqWeightExponent
        namelist /specPartition/ shapeAdaptationEnabled, scaleOptimizationEnabled, shapeOptimizationEnabled

        assertion = .true._LK

        rinitEnabled = compiler_version() == "GCC version 10.2.0"

        ! read the number of tests

        ntest = 1
        open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
        read(test%File%unit, nml = specTest)
        close(test%File%unit)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do itest = 1, ntest

            ! read cluster spec

            relTol = 0.01_RK
            isRepeatable = .false._LK
            rngseed = -huge(rngseed)
            if (.not. allocated(dist)) allocate(character(63) :: dist)
            open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
            read(test%File%unit, nml = specData)
            dist = trim(adjustl(dist))
            close(test%File%unit)
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true._LK, image_distinct = .true._LK)
            else
                if(rngseed /= -huge(rngseed)) randomSeed = randomSeed_type(imageID = test%image%id, inputSeed = rngseed)
            end if
!nc = 4
!sizeMin = 1000
!sizeMax = 1000
!centerMin = -5._RK
!centerMax = 5._RK
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
                                    , hubEnabled = .false._LK & ! LCOV_EXCL_LINE
                                    , relTol = relTol & ! LCOV_EXCL_LINE
                                    )

            ! write data to output for further investigation

            test%File = test%openFile(label = "ClusteredPoint")
            call ClusteredPoint%write(test%File%unit)
            close(test%File%unit)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Partition the clustered points clustered points
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! read partition spec

            nsim = 0
            minSize = 0 !nd + 1
            stanEnabled = .true._LK
            isRepeatable = .false._LK
            rngseed = -huge(rngseed)
            expansionMinVol = 1._RK
            expansionMaxDen = 1._RK
            shrinkageMinVol = 1._RK
            shrinkageMaxDen = 1._RK
            mahalSqWeightExponent = 0._RK
            shapeAdaptationEnabled = .true._LK
            scaleOptimizationEnabled = .true._LK
            shapeOptimizationEnabled = .true._LK
            kvolumeNumRecursionMax = 100
            nemax = ClusteredPoint%np / (ClusteredPoint%nd + 1)
            open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_partition@test_constructPartition_1.nml", status = "old")
            read(test%File%unit, nml = specPartition)
            close(test%File%unit)
            !if (shrinkageMaxDen > 0._RK) then; shrinkageMaxDen = log(shrinkageMaxDen); else; shrinkageMaxDen = NEGBIG_RK; endif ! LCOV_EXCL_LINE
            if (isRepeatable .and. rinitEnabled) then
                call random_init (repeatable = .true._LK, image_distinct = .true._LK)
            else
                if(rngseed /= -huge(rngseed)) randomSeed = randomSeed_type(imageID = test%image%id, inputSeed = rngseed)
            end if

#if BENCHM_ENABLED

            Partition = getPartition( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                    , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! & ! LCOV_EXCL_LINE
                                    , minSize = minSize & ! LCOV_EXCL_LINE
                                    , nemax = nemax & ! LCOV_EXCL_LINE
                                    )

#else

            !write(*,*) "zeroth nemax", nemax, ClusteredPoint%np, ClusteredPoint%nd + 1
            Partition = partition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                        , nc = nc & ! LCOV_EXCL_LINE
                                        , nt = nt & ! LCOV_EXCL_LINE
                                        , nsim = nsim & ! LCOV_EXCL_LINE
                                        , nemax = nemax & ! LCOV_EXCL_LINE
                                        , trimEnabled = .false._LK & ! LCOV_EXCL_LINE
!#if defined MINVOL
!                                        , stanEnabled = .false._LK & ! LCOV_EXCL_LINE
!#endif
                                        , stanEnabled = stanEnabled & ! LCOV_EXCL_LINE
#if defined MINVOL_RECURSIVE
                                        , logExpansion = log(expansionMinVol) & ! LCOV_EXCL_LINE
                                        , logShrinkage = log(shrinkageMinVol) & ! LCOV_EXCL_LINE
                                        , kvolumeNumRecursionMax = min(kvolumeNumRecursionMax, 20) & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = min(mahalSqWeightExponent,1._RK) & ! LCOV_EXCL_LINE
#elif defined MAXDEN_RECURSIVE || defined MAXDEN
                                        , logExpansion = log(expansionMaxDen) & ! LCOV_EXCL_LINE
                                        , logShrinkage = log(shrinkageMaxDen) & ! LCOV_EXCL_LINE
                                        , scaleOptimizationEnabled = scaleOptimizationEnabled & ! LCOV_EXCL_LINE
                                        , shapeOptimizationEnabled = shapeOptimizationEnabled & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#if defined MAXDEN_RECURSIVE
                                        , kvolumeNumRecursionMax = min(kvolumeNumRecursionMax, 20) & ! LCOV_EXCL_LINE
                                        !, kvolumeNumRecursionMax = kvolumeNumRecursionMax & ! LCOV_EXCL_LINE
#else
                                        , kvolumeNumRecursionMax = kvolumeNumRecursionMax & ! LCOV_EXCL_LINE
                                        , shapeAdaptationEnabled = shapeAdaptationEnabled & ! LCOV_EXCL_LINE
#endif
#endif
#if defined MAXDEN_RECURSIVE || defined MINVOL_RECURSIVE
                                        , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
#endif
                                        , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! & ! LCOV_EXCL_LINE
                                       !, nemax = ClusteredPoint%nemax & ! LCOV_EXCL_LINE
                                       !, kmeansNumFailMax = 10000 & ! LCOV_EXCL_LINE
                                       !, kmeansNumRecursionMax = 10000 & ! LCOV_EXCL_LINE
                                       !, mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
                                        , minSize = minSize & ! LCOV_EXCL_LINE
                                        )

            ! write data to output for further investigation

            test%File = test%openFile(label = "Partition")
            call partition%write(test%File%unit, ClusteredPoint%Point)
            close(test%File%unit)

            assertion = assertion .and. .not. partition%err%occurred
            call test%assert(assertion)
            assertion = assertion .and. partition%err%stat /= 1
            call test%assert(assertion)
            assertion = assertion .and. partition%err%stat /= 2
            call test%assert(assertion)
            assertion = assertion .and. all(partition%Size(1:partition%ne) >= 0)
            call test%assert(assertion)
            assertion = assertion .and. all(partition%membership > 0) .and. all(partition%membership < partition%ne + 1)
            call test%assert(assertion)

            !if (.not. assertion) return

            ! Ensure all points are within their corresponding clusters.

            do ip = 1, partition%np
                ic = partition%membership(ip)
                NormedPoint = ClusteredPoint%Point(:,ip) - partition%Center(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(partition%invCov(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    if (test%traceable) write(test%disp%unit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                    assertion = assertion .and. isInside
                    exit
                end if
            end do

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne                 =", partition%ne
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%Size < 1           =", pack(partition%Size(1:partition%ne), mask = partition%Size(1:partition%ne) < 1)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership < 1     =", pack(partition%membership, mask = partition%membership < 1)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%membership > ne    =", pack(partition%membership, mask = partition%membership > partition%ne)
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%occurred       =", partition%err%occurred
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%stat           =", partition%err%stat
                write(test%disp%unit,"(*(g0.15,:,' '))") "partition%err%msg            =", partition%err%msg
                if (.not. isInside) write(test%disp%unit,"(*(g0.15,:,' '))") "One or more points are NOT bounded!"
                write(test%disp%unit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
!assertion = .true._LK
            call test%assert(assertion)

            ! The membership IDs must cover a full range from 1 to partition%ne

            block

                use pm_arrayUnique, only: setUnique
                integer(IK) :: lenUnique
                integer(IK), allocatable :: UniqueValue(:), UniqueCount(:)

                call setUnique  ( Array = partition%membership & ! LCOV_EXCL_LINE
                                , Unique = UniqueValue & ! LCOV_EXCL_LINE
                                , Count = UniqueCount & ! LCOV_EXCL_LINE
                                )
                lenUnique = size(UniqueValue, kind = IK)
                do ic = 1, partition%ne
                    assertion = assertion .and. any(UniqueValue == ic)
                    if (test%traceable .and. .not. assertion) then
                        ! LCOV_EXCL_START
                        write(test%disp%unit,"(*(g0.15,:,' '))")
                        write(test%disp%unit,"(*(g0.15,:,' '))") "Membership IDs must cover a full range 1:partition%ne"
                        write(test%disp%unit,"(*(g0.15,:,' '))") "partition%ne             =", partition%ne
                        write(test%disp%unit,"(*(g0.15,:,' '))") "Unique Membership IDs    =", UniqueValue
                        write(test%disp%unit,"(*(g0.15,:,' '))")
                        ! LCOV_EXCL_STOP
                    end if
                    if (.not. assertion) return
                end do

            end block

            if (.not. assertion) return
#endif

        end do

    end function test_constructPartition_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MAXDEN
    !> Output the bias predictions for comparison with MATLAB.
    function test_getVarCorrectionScaleFactorSq_1() result(assertion)
        implicit none
        logical(LK)             :: assertion
        integer(IK)             :: ip, id, npmin, npmax
        real(RK), allocatable   :: BiasCorrectionScaleFactorSq(:)
        assertion = .true._LK
        test%File = test%openFile(label = "BiasCorrection")
        do id = 2, 31
            npmin = id + 1_IK
            npmax = 20000_IK + id
            if (allocated(BiasCorrectionScaleFactorSq)) deallocate(BiasCorrectionScaleFactorSq)
            allocate(BiasCorrectionScaleFactorSq(npmin:npmax))
            BiasCorrectionScaleFactorSq(npmin:npmax) = getVarCorrectionScaleFactorSq(id, npmin, npmax)
            write(test%File%unit, "(*(g0,:,','))") id, sqrt(BiasCorrectionScaleFactorSq)
        end do
        close(test%File%unit)
    end function test_getVarCorrectionScaleFactorSq_1
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
