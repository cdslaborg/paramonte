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

        use iso_fortran_env, only: compiler_version
        use ClusteredPoint_pmod, only: ClusteredPoint_type
        use pm_randomSeed, only: randomSeed_type
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
#if MAXDEN_ENABLED  
        type(PartitionMaxDen_type)      :: Partition
#elif MINVOL_ENABLED    
        type(PartitionMinVol_type)      :: Partition
#endif  
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

        integer(IK)                     :: nt, nemax, minSize
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
        namelist /specPartition/ rngseed, isRepeatable, nc, nt, nemax, minSize, inclusionFraction, stanEnabled
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

            minSize = 0 !nd + 1
            stanEnabled = .true._LK
            isRepeatable = .false._LK
            rngseed = -huge(rngseed)
            expansionMinVol = 2._RK
            expansionMaxDen = 1._RK
            shrinkageMinVol = 2._RK
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

#if MAXDEN_ENABLED
            Partition = PartitionMaxDen_type( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
#elif MINVOL_ENABLED
            Partition = PartitionMinVol_type( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
#endif
                                            , nc = nc & ! LCOV_EXCL_LINE
                                            , nt = nt & ! LCOV_EXCL_LINE
                                            , nemax = nemax & ! LCOV_EXCL_LINE
                                            , trimEnabled = .false._LK & ! LCOV_EXCL_LINE
                                            , stanEnabled = stanEnabled & ! LCOV_EXCL_LINE
#if MAXDEN_ENABLED
                                            , scaleOptimizationEnabled = scaleOptimizationEnabled & ! LCOV_EXCL_LINE
                                            , shapeOptimizationEnabled = shapeOptimizationEnabled & ! LCOV_EXCL_LINE
                                            , shapeAdaptationEnabled = shapeAdaptationEnabled & ! LCOV_EXCL_LINE
                                            , logExpansion = log(expansionMaxDen) & ! LCOV_EXCL_LINE
                                            , logShrinkage = log(shrinkageMaxDen) & ! LCOV_EXCL_LINE
#elif MINVOL_ENABLED
                                            , logExpansion = log(expansionMinVol) & ! LCOV_EXCL_LINE
                                            , logShrinkage = log(shrinkageMinVol) & ! LCOV_EXCL_LINE
#if DYNESTY_ENABLED
                                            , method = "dynesty" & ! LCOV_EXCL_LINE
#elif MULTINEST_ENABLED
                                            , method = "multinest" & ! LCOV_EXCL_LINE
#endif
#endif
                                            , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                            , kvolumeNumRecursionMax = kvolumeNumRecursionMax & ! LCOV_EXCL_LINE
                                            , biasCorrectionEnabled = biasCorrectionEnabled & ! LCOV_EXCL_LINE
                                            , parentLogVolNormed = ClusteredPoint%logSumVolNormedEffective & ! & ! LCOV_EXCL_LINE
                                            , minSize = minSize & ! LCOV_EXCL_LINE
                                            )

            ! write data to output for further investigation

            test%File = test%openFile(label = "Partition")
            call partition%write(test%File%unit, ClusteredPoint%Point)
            close(test%File%unit)

            assertion = assertion .and. .not. partition%err%occurred
!write(*,*) assertion
            assertion = assertion .and. partition%err%stat /= 1
!write(*,*) assertion
            assertion = assertion .and. partition%err%stat /= 2
!write(*,*) assertion
            assertion = assertion .and. all(partition%Size(1:partition%ne) >= 0)
!write(*,*) assertion
            assertion = assertion .and. all(partition%membership > 0) .and. all(partition%membership < partition%ne + 1)
!write(*,*) assertion

            !if (.not. assertion) return

            ! Ensure all points are within their corresponding clusters.

            do ip = 1, partition%np
                ic = partition%membership(ip)
                NormedPoint = ClusteredPoint%Point(:,ip) - partition%Center(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(partition%invCov(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    if (test%traceable) then
                        write(test%disp%unit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                        write(*,*) assertion, mahalSq
                    end if
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

        end do