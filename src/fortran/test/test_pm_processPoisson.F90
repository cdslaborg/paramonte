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

!>  \brief This module contains tests of the module [pm_procPois](@ref pm_procPois).
!>  \author Amir Shahmoradi

module test_pm_processPoisson

    use pm_procPois
    use pm_test, only: test_type, LK
    use pm_container, only: IV => cvi_type, RV => cvr_type

    implicit none

    private
    public :: setTest
    type(test_type) :: test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 10
        integer(IK)             :: HubMinDistEdgeIndex(10) =    [5, 8, 7, 7, 10, 9, 3, 2, 3, 5]
        real(RK)                :: HubMinDistEdgeLenSq(10) =    [ 0.693595837531895E-1_RK &
                                                                , 0.411631653715387E-2_RK &
                                                                , 0.324022691025059E-2_RK &
                                                                , 0.122396412769469E-1_RK &
                                                                , 0.106512144573307E-1_RK &
                                                                , 0.442485474967571000_RK &
                                                                , 0.324022691025059E-2_RK &
                                                                , 0.411631653715387E-2_RK &
                                                                , 0.375960647686508E-1_RK &
                                                                , 0.106512144573307E-1_RK &
                                                                ]
        integer(IK)             :: HubNodeIndex(7) = [3, 7, 5, 2, 9, 10, 8]
        integer(IK)             :: HubEdgeCount(7) = [2, 2, 2, 1, 1, 1, 1]
        type(IV), allocatable   :: HubEdgeIndex(:)
        type(RV), allocatable   :: HubEdgeLenSq(:)
        real(RK)                :: Point(2,10) = reshape( [ 0.278498218867048_RK, 0.421761282626275_RK &
                                                        , 0.546881519204984_RK, 0.915735525189067_RK &
                                                        , 0.957506835434298_RK, 0.792207329559554_RK &
                                                        , 0.964888535199277_RK, 0.959492426392903_RK &
                                                        , 0.157613081677548_RK, 0.655740699156587_RK &
                                                        , 0.970592781760616_RK, 0.035711678574190_RK &
                                                        , 0.957166948242946_RK, 0.849129305868777_RK &
                                                        , 0.485375648722841_RK, 0.933993247757551_RK &
                                                        , 0.800280468888800_RK, 0.678735154857773_RK &
                                                        , 0.141886338627215_RK, 0.757740130578333_RK &
                                                        ], shape = [2, 10] )
        real(RK)                :: Dist(10,10) = reshape( [ 0._RK &
                                                        , 0.562174482003378_RK &
                                                        , 0.773487540339897_RK &
                                                        , 0.871944063189390_RK &
                                                        , 0.263362077287504_RK &
                                                        , 0.792482921440967_RK &
                                                        , 0.802019121669115_RK &
                                                        , 0.552430861815293_RK &
                                                        , 0.581628994675654_RK &
                                                        , 0.362690766485521_RK &
                                                        , 0.562174482003378_RK &
                                                        , 0._RK &
                                                        , 0.428803411185017_RK &
                                                        , 0.420291008496988_RK &
                                                        , 0.468110271216848_RK &
                                                        , 0.976715518780844_RK &
                                                        , 0.415656735459690_RK &
                                                        , 0.064158526613022_RK &
                                                        , 0.346958503625479_RK &
                                                        , 0.434722487351898_RK &
                                                        , 0.773487540339897_RK &
                                                        , 0.428803411185017_RK &
                                                        , 0._RK &
                                                        , 0.167447881784044_RK &
                                                        , 0.811451266874729_RK &
                                                        , 0.756608823601091_RK &
                                                        , 0.056922991051512_RK &
                                                        , 0.492961564490394_RK &
                                                        , 0.193897046828080_RK &
                                                        , 0.816348444365176_RK &
                                                        , 0.871944063189390_RK &
                                                        , 0.420291008496988_RK &
                                                        , 0.167447881784044_RK &
                                                        , 0._RK &
                                                        , 0.862530445641055_RK &
                                                        , 0.923798359204721_RK &
                                                        , 0.110632912268217_RK &
                                                        , 0.480190395997296_RK &
                                                        , 0.325454237972598_RK &
                                                        , 0.847370405683894_RK &
                                                        , 0.263362077287504_RK &
                                                        , 0.468110271216848_RK &
                                                        , 0.811451266874729_RK &
                                                        , 0.862530445641055_RK &
                                                        , 0._RK &
                                                        , 1.022434339755630_RK &
                                                        , 0.822608982898776_RK &
                                                        , 0.429945090865161_RK &
                                                        , 0.643078623169773_RK &
                                                        , 0.103204721100010_RK &
                                                        , 0.792482921440967_RK &
                                                        , 0.976715518780844_RK &
                                                        , 0.756608823601091_RK &
                                                        , 0.923798359204721_RK &
                                                        , 1.02243433975563_RK &
                                                        , 0._RK &
                                                        , 0.813528419539970_RK &
                                                        , 1.020953203495600_RK &
                                                        , 0.665195817009978_RK &
                                                        , 1.099126678046850_RK &
                                                        , 0.802019121669115_RK &
                                                        , 0.415656735459690_RK &
                                                        , 0.056922991051512_RK &
                                                        , 0.110632912268217_RK &
                                                        , 0.822608982898776_RK &
                                                        , 0.813528419539970_RK &
                                                        , 0._RK &
                                                        , 0.479363034594627_RK &
                                                        , 0.231619373332412_RK &
                                                        , 0.820386770843889_RK &
                                                        , 0.552430861815293_RK &
                                                        , 0.064158526613022_RK &
                                                        , 0.492961564490394_RK &
                                                        , 0.480190395997296_RK &
                                                        , 0.429945090865161_RK &
                                                        , 1.020953203495600_RK &
                                                        , 0.479363034594627_RK &
                                                        , 0._RK &
                                                        , 0.405366179835696_RK &
                                                        , 0.386070029224440_RK &
                                                        , 0.581628994675654_RK &
                                                        , 0.346958503625479_RK &
                                                        , 0.193897046828080_RK &
                                                        , 0.325454237972598_RK &
                                                        , 0.643078623169773_RK &
                                                        , 0.665195817009978_RK &
                                                        , 0.231619373332412_RK &
                                                        , 0.405366179835696_RK &
                                                        , 0._RK &
                                                        , 0.663117347798649_RK &
                                                        , 0.362690766485521_RK &
                                                        , 0.434722487351898_RK &
                                                        , 0.816348444365176_RK &
                                                        , 0.847370405683894_RK &
                                                        , 0.103204721100010_RK &
                                                        , 1.099126678046850_RK &
                                                        , 0.820386770843889_RK &
                                                        , 0.386070029224440_RK &
                                                        , 0.663117347798649_RK &
                                                        , 0._RK &
                                                        ], shape = [10,10] )
    end type TestData_type

   !type(TestData_type) :: TestData

    interface TestData_type
        module procedure :: contructTestData
    end interface TestData_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        !TestData = TestData_type()
        call test%run(test_getLogDensity_1, SK_"test_getLogDensity_1")
        call test%run(test_getPosteriorDistSortedExpDiff_1, SK_"test_getPosteriorDistSortedExpDiff_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function contructTestData() result(TestData)
        implicit none
        type(TestData_type) :: TestData
        allocate(TestData%HubEdgeIndex(7))
        TestData%HubEdgeIndex(1)%val = [7, 9]
        TestData%HubEdgeIndex(2)%val = [3, 4]
        TestData%HubEdgeIndex(3)%val = [1,10]
        TestData%HubEdgeIndex(4)%val = [8]
        TestData%HubEdgeIndex(5)%val = [6]
        TestData%HubEdgeIndex(6)%val = [5]
        TestData%HubEdgeIndex(7)%val = [2]
        allocate(TestData%HubEdgeLenSq(7))
        TestData%HubEdgeLenSq(1)%val = [0.324022691025059E-2_RK, 0.375960647686508E-1_RK]
        TestData%HubEdgeLenSq(2)%val = [0.324022691025059E-2_RK, 0.122396412769469E-1_RK]
        TestData%HubEdgeLenSq(3)%val = [0.693595837531895E-1_RK, 0.106512144573307E-1_RK]
        TestData%HubEdgeLenSq(4)%val = [0.411631653715387E-2_RK]
        TestData%HubEdgeLenSq(5)%val = [0.442485474967571_RK]
        TestData%HubEdgeLenSq(6)%val = [0.106512144573307E-1_RK]
        TestData%HubEdgeLenSq(7)%val = [0.411631653715387E-2_RK]
    end function contructTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getLogExpectedMinPairDist()`.
    function test_getLogDensity_1() result(assertion)

        use ClusteredPoint_pmod, only: ClusteredPoint_type
        use pm_ellipsoid, only: getLogVolUnitBall
        use pm_kind, only: IK, RK

        implicit none
        logical(LK)                     :: assertion
    
        type(ClusteredPoint_type)       :: ClusteredPoint
        real(RK)                        :: logDensityFromExpectedMinPairDist
        real(RK)                        :: logVolUnitBall
        real(RK)                        :: logExpectedMinPairDist
        integer(IK)                     :: i

        character(:, SK), allocatable   :: dist
        integer(IK)                     :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                        :: etamin, etamax, centerMin, centerMax
        namelist /specData/ rngseed, nd, nc, sizeMin, sizeMax, etamin, etamax, centerMin, centerMax, dist

        assertion = .true._LK

        rngseed = -huge(rngseed)
        allocate(character(63) :: dist)
        open(newunit = test%File%unit, file = test%dir%inp//"/test_pm_processPoisson@test_getLogDensity_1.nml", status = "old")
        read(test%File%unit, nml = specData)
        dist = trim(adjustl(dist))
        close(test%File%unit)

        logVolUnitBall = getLogVolUnitBall(nd)

!write(*,*) "logVolUnitBall", logVolUnitBall

        test%File = test%openFile() ! label = "LogVolEstimate")

        do i = 1, 5 ! 000

!if (mod(i,100)==0) write(*,*) "i = ", i

            ! Generate clustered points

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
                                    , hubEnabled = .true. & ! LCOV_EXCL_LINE
                                    )

            ! write data to output for further investigation

            !test%File = test%openFile(label = "ClusteredPoint")
            !call ClusteredPoint%write(test%File%unit)
            !close(test%File%unit)

            !logExpectedMinPairDist = 0._RK
            !nh = ClusteredPoint%Hub%nh / 2
            !do ih = 1, nh
            !    logExpectedMinPairDist = logExpectedMinPairDist + sum(sqrt(ClusteredPoint%Hub%EdgeLenSq(ih)%val(1:1)))
            !    if (size(ClusteredPoint%Hub%EdgeLenSq(ih)%val) /= ClusteredPoint%Hub%EdgeCount(ih)) error stop
            !end do
            !logExpectedMinPairDist = log( logExpectedMinPairDist / nh ) ! / ClusteredPoint%Hub%EdgeCount(1:2)

            !logExpectedMinPairDist = log( sum(sqrt(ClusteredPoint%Hub%EdgeLenSq(1)%val)) / ClusteredPoint%Hub%EdgeCount(1) )
            !if (ClusteredPoint%Hub%EdgeCount(1) /= size((ClusteredPoint%Hub%EdgeLenSq(1)%val))) error stop

            logExpectedMinPairDist = log( sqrt(minval(ClusteredPoint%Hub%EdgeLenSq(1)%val)) )

            logDensityFromExpectedMinPairDist = getLogDensity(nd = nd, logVolUnitBall = logVolUnitBall, logExpectedMinPairDist = logExpectedMinPairDist)

            write(test%File%unit, "(*(g0,:,','))") ClusteredPoint%LogDensity, logDensityFromExpectedMinPairDist

        end do

        close(test%File%unit)

    end function test_getLogDensity_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getPosteriorDistSortedExpDiff()` with a uniform cubic distribution.
    function test_getPosteriorDistSortedExpDiff_1() result(assertion)

        use pm_arraySpace, only: getLogSpace, getLinSpace
        use pm_knn, only: setDisSortedExpDiff
        use pm_ellipsoid, only: getLogVolUnitBall
        use pm_domainBall, only: getUnifRand
        use pm_distExp, only: setExpRand
        use pm_kind, only: IK, RK
        use pm_arrayUnique, only: getUnique
        use pm_err, only: err_type

        use pm_sampleWeight, only: setChoLowCovUpp
        use pm_distanceMahal, only: getDisMahalSq
        use pm_matrixDet, only: getMatInvFromChoLow
        use pm_matrixChol, only: setChoLow
        use pm_matrixInit, only: setMatInit

        logical(LK)                 :: assertion
        integer(IK) , parameter     :: nd = 5_IK, np = 10_IK**4 !* 2
        real(RK)                    :: Sample(nd,np)
        real(RK)                    :: Reference(nd)
        real(RK)                    :: disSortedExpDiff(np)
        real(RK)                    :: ProbKS(np),StatKS(np),Lambda(np)
        integer(IK)                 :: ShuffleIndex(np)
        integer(IK), allocatable    :: SortedIndex(:)
        real(RK)                    :: logVolUnitBall
        real(RK)                    :: alpha, beta
        real(RK)                    :: volumeAvg, volumeStd
        real(RK)                    :: densityAvg, densityStd
        logical(LK)                 :: failed
        integer(IK)                 :: divisor
       !integer(IK)                 :: ip,i,j
        type(err_type)              :: Err

        real(RK)                    :: choLowCovUpp(nd,nd)
        real(RK)                    :: choDia(nd)
        real(RK)                    :: mean(nd)
        real(RK)                    :: volumeRef
        real(RK)                    :: logVol

        call random_seed()
        !call random_init(repeatable = .true._LK, image_distinct = .true._LK)
        !block
        !    use pm_randomSeed, only: randomSeed_type
        !    type(randomSeed_type) :: randomSeed
        !    randomSeed = randomSeed_type(imageID = 1_IK, inputSeed = 1234_IK)
        !end block

        !block
        !use pm_arraySpace, only: getLinSpace
        !write(*,*) getLinSpace(1._RK, real(10,RK), 10)
        !end block

        assertion = .true._LK

        logVolUnitBall = getLogVolUnitBall(nd)

        ! Generate points uniformly distributed within a cube.

        !call random_number(Sample)
        !Reference(1:nd) = 0.5_RK
        Reference(1:nd) = -10.0_RK

        ! Generate points uniformly distributed within an ellipsoid.

       !choLowCovUpp = getVecDia(nd, nd, [(real(i,RK),i=1,nd)], 0.8_RK)
        call setMatInit(choLowCovUpp, 1._RK)
        call setChoLow(choLowCovUpp, choDia, nd)
        assertion = choDia(1) > 0._RK
        !write(*,*) choDia(1)
        !write(*,*) ((choLowCovUpp(i,j), j = 1, nd), new_line("a"), i = 1, nd)
        call test%assert(assertion)
        !choDia = exp((-logVolUnitBall)/nd) ! - sum(log(choDia))

        divisor = 2 ! This must be set such that the different sections of the volume have the same uniform distribution.
        mean = -10._RK
        call getUnifRand(rand = Sample(1:nd,1:np/divisor), mean = mean, choLow = choLowCovUpp, choDia = choDia)
        mean = +10._RK
        call getUnifRand(rand = Sample(1:nd,np/divisor+1:np), mean = mean, choLow = choLowCovUpp, choDia = choDia)

        logVol = sum(log(choDia)) + logVolUnitBall + log(2._RK)
        volumeRef = exp(logVol)
        write(*, "(*(g0.8,:,', '))") "choDia: ", choDia
        write(*, "(*(g0.8,:,', '))") "logVolUnitBall", logVolUnitBall
        write(*, "(*(g0.8,:,', '))") "volUnitBall", exp(logVolUnitBall)
        write(*, "(*(g0.8,:,', '))") "volume: ", volumeRef

#define MAHAL_TRANSFORM_ENABLED 1

#if     MAHAL_TRANSFORM_ENABLED
        block
            use pm_sampleWeight, only: transform
            use pm_matrixDet, only: getLogPDF
            real(RK)    :: InvMatLow(nd,nd)!, invCov(nd,nd), mahalSq(np)

            mean = sum(Sample, dim = 2) / np
            call setChoLowCovUpp(nd,np,mean,Sample,choLowCovUpp,choDia)
            assertion = choDia(1) > 0._RK
            call test%assert(assertion)

            ! Compute the correlation matrix.

            InvMatLow = getLogPDF(nd,choLowCovUpp,choDia)
            !invCov = getMatInvFromChoLow(nd,choLowCovUpp,choDia)
            !mahalSq = getDisMahalSq(nd,np,mean,invCov,Sample)

            ! Decorrelate the entire sample.

            Sample(1:nd,1:np) = transform(Sample = Sample, choLow = InvMatLow)
        end block
#endif

        ! Compute the differences of the distances of points from each other.

        call setDisSortedExpDiff(Sample, Reference, disSortedExpDiff, ShuffleIndex)
        if (failed) then; assertion = .false.; return; end if ! LCOV_EXCL_LINE

        ! Write points to output file.

        !block
        !    integer(IK) :: fileUnit
        !    open(newunit = fileUnit, file = "Points.txt", status = "replace")
        !    write(fileUnit,*) (Sample(1:nd,ip), new_line("a"), ip = 1, np)
        !    close(fileUnit)
        !end block

        ! Get true random exponentials.

        !disSortedExpDiff = [(genExpRand(mean = exp(-logVolUnitBall)), ip = 1, np)]

        !SortedIndex = [(ip, ip = 1, np)]
        !SortedIndex = [(ip, ip = 1, np, np/100)]
        !SortedIndex = getUnique(nint(getLinSpace(1._RK, real(np,RK), np), IK))
        SortedIndex = getUnique(nint(getLogSpace(0._RK, log(real(np,RK)), 10000_IK), kind = IK))
        !write(*,*) SortedIndex
        call getDistSortedExpDiffWeightKS2  ( SortedIndex = SortedIndex & ! LCOV_EXCL_LINE
                                            , disSortedExpDiff = disSortedExpDiff & ! LCOV_EXCL_LINE
                                            , logVolUnitBall = logVolUnitBall & ! LCOV_EXCL_LINE
                                            , ProbKS = ProbKS & ! LCOV_EXCL_LINE
                                            , StatKS = StatKS & ! LCOV_EXCL_LINE
                                            , Lambda = Lambda & ! LCOV_EXCL_LINE
                                            , alpha = alpha & ! LCOV_EXCL_LINE
                                            , beta = beta & ! LCOV_EXCL_LINE
                                            )
        if (err%occurred) then; assertion = .false.; return; end if ! LCOV_EXCL_LINE

#if     MAHAL_TRANSFORM_ENABLED
        beta = exp(log(beta) - sum(log(choDia)))
#endif
        densityAvg = alpha / beta
        densityStd = sqrt(alpha) / beta
        volumeAvg = np * beta / (alpha - 1._RK) ! np / densityAvg
        volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! np / densityStd
        write(*, "(*(g0.8,:,', '))")
        write(*, "(*(g0.8,:,', '))") "Weighted Likelihood"
        write(*, "(*(g0.8,:,', '))") "nd, np", nd, np
        write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
        write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
        write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
        write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd

        block
            integer(IK) :: effSamSize
            effSamSize = nint(alpha,IK)
            alpha = real(effSamSize,RK)
            beta = exp(log(sum(disSortedExpDiff(1:effSamSize))) + logVolUnitBall)
#if         MAHAL_TRANSFORM_ENABLED
            beta = exp(log(beta) - sum(log(choDia)))
#endif
            densityAvg = alpha / beta
            densityStd = sqrt(alpha) / beta
            volumeAvg = np * beta / (alpha - 1._RK) ! np / densityAvg
            volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! np / densityStd
            write(*, "(*(g0.8,:,', '))")
            write(*, "(*(g0.8,:,', '))") "sum(Weight) as effective sample size"
            write(*, "(*(g0.8,:,', '))") "nd, np", nd, np
            write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
            write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd
        end block

#if     CHECK_ENABLED

        ! Write the weights to output files.

        block
            use iso_fortran_env, only: output_unit
            use pm_arraySort, only: setSorted
            use pm_mathCumSum, only: getCumSum
            integer(IK)                 :: ip
            integer(IK)                 :: fileUnit
            real(RK)                    :: BetaArray(np)
            real(RK)                    :: CumSumDistSortedExpDiff(np)
            real(RK)                    :: WeightedDistSortedExpDiff(np)
            CumSumDistSortedExpDiff = getCumSum(disSortedExpDiff)
            WeightedDistSortedExpDiff = ProbKS * disSortedExpDiff
            BetaArray = exp(log(CumSumDistSortedExpDiff) + logVolUnitBall - logVol)
#if         MAHAL_TRANSFORM_ENABLED
            BetaArray = BetaArray / product(choDia)
#endif

            open(newunit = fileUnit, file = "WeightsAll.txt", status = "replace")
            write(fileUnit,"(*(g0,:,','))") "ProbKS, StatKS, Lambda, disSortedExpDiff, WeightedDistSortedExpDiff, VolumeAvg, VolumeStd"
            write(fileUnit,"(*(g0,:,','))") 0._RK, 0._RK, 0._RK, 0._RK, 0._RK, exp(logVol-logVol), 0._RK
            write(fileUnit,"(*(g0))")(ProbKS(ip), "," & ! LCOV_EXCL_LINE
                                    , StatKS(ip), "," & ! LCOV_EXCL_LINE
                                    , Lambda(ip), "," & ! LCOV_EXCL_LINE
                                    , disSortedExpDiff(ip), "," & ! LCOV_EXCL_LINE
                                    , WeightedDistSortedExpDiff(ip), "," & ! LCOV_EXCL_LINE
                                    , np * BetaArray(ip) / (ip - 1._RK), "," & ! LCOV_EXCL_LINE
                                    , np * BetaArray(ip) / (ip - 1._RK) / sqrt(ip - 2._RK) & ! LCOV_EXCL_LINE ! volumeStd
                                    , new_line("a"), ip = 2, size(disSortedExpDiff))
            write(*,*) "maxloc(ProbKS) = ", maxloc(ProbKS(2:) + 1)
        end block

        ! Identify the last point where its likelihood of occurrence is above the expected value.

        block
            integer(IK) :: effSamSize
            real(RK) :: ExpectedSampleSize(np)
            where(ProbKS > 0._RK)
                ExpectedSampleSize = 1._RK / ProbKS
            elsewhere
                ExpectedSampleSize = huge(1._RK)
            endwhere
            do effSamSize = np, 1, -1
                if (ExpectedSampleSize(effSamSize) <= np - effSamSize + 1) exit
            end do

            alpha = effSamSize
            beta = exp(log(sum(disSortedExpDiff(1:effSamSize))) + logVolUnitBall)
#if MAHAL_TRANSFORM_ENABLED
            beta = exp(log(beta) - sum(log(choDia)))
#endif
            densityAvg = alpha / beta
            densityStd = sqrt(alpha) / beta
            volumeAvg = np * beta / (alpha - 1._RK) ! effSamSize / densityAvg
            volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! effSamSize / densityStd
            write(*, "(*(g0.8,:,', '))")
            write(*, "(*(g0.8,:,', '))") "ProbKS > 1 / SampleSize"
            write(*, "(*(g0.8,:,', '))") "effSamSize", effSamSize
            write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
            write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd

            effSamSize = int(0.5_RK * effSamSize, IK)
            alpha = effSamSize
            beta = exp(log(sum(disSortedExpDiff(1:effSamSize))) + logVolUnitBall)
#if MAHAL_TRANSFORM_ENABLED
            beta = exp(log(beta) - sum(log(choDia)))
#endif
            densityAvg = alpha / beta
            densityStd = sqrt(alpha) / beta
            volumeAvg = np * beta / (alpha - 1._RK) ! effSamSize / densityAvg
            volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! effSamSize / densityStd
            write(*, "(*(g0.8,:,', '))")
            write(*, "(*(g0.8,:,', '))") "half ProbKS > 1 / SampleSize"
            write(*, "(*(g0.8,:,', '))") "effSamSize", effSamSize
            write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
            write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd

            !effSamSize = 300_IK
            !alpha = effSamSize
            !beta = exp( log(sum(disSortedExpDiff(1:effSamSize))) + logVolUnitBall )
            !densityAvg = alpha / beta
            !densityStd = sqrt(alpha) / beta
            !volumeAvg = np * beta / (alpha - 1._RK) ! effSamSize / densityAvg
            !volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! effSamSize / densityStd
            !write(*, "(*(g0.8,:,', '))")
            !write(*, "(*(g0.8,:,', '))") "effSamSize", effSamSize
            !write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            !write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            !write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
        end block

        block
            use pm_mathCumSum, only: getCumSum
            real(RK), allocatable :: CumSumReverseAlpha(:)
           !CumSumReverseAlpha = getCumSumReverse(np,ProbKS)
            CumSumReverseAlpha = getCumSum(ProbKS, backward = .true._LK, reversed = .true._LK)
            CumSumReverseAlpha = np * CumSumReverseAlpha / sum(CumSumReverseAlpha)
            beta = exp(log(sum(disSortedExpDiff * CumSumReverseAlpha)) + logVolUnitBall)
#if MAHAL_TRANSFORM_ENABLED
            beta = exp(log(beta) - sum(log(choDia)))
#endif
            alpha = sum(CumSumReverseAlpha)
            densityAvg = alpha / beta
            densityStd = sqrt(alpha) / beta
            volumeAvg = np * beta / (alpha - 1._RK) ! effSamSize / densityAvg
            volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! effSamSize / densityStd
            write(*, "(*(g0.8,:,', '))")
            write(*, "(*(g0.8,:,', '))") "Weighted Samples"
            write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
            write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd
        end block

        block
            use pm_arraySpace, only: getLinSpace
            integer(IK) :: minLocStatKS
            !real(RK), allocatable :: CorrectedStatKS(:)
            !!minLocStatKS = minloc(StatKS * sqrt([(real(ip,RK), ip = 1, np)]), dim = 1)
            !CorrectedStatKS = StatKS * sqrt(getLinSpace(1._RK, real(np,RK), np))
            !write(*,*) real(getLinSpace(1._RK, real(100,RK), 100))
            !write(*,*) real(CorrectedStatKS(1:100))
            !minLocStatKS = minloc(CorrectedStatKS(2:), dim = 1) + 1
            minLocStatKS = minloc(Lambda(2:), dim = 1) + 1
            beta = exp(log(sum(disSortedExpDiff(1:minLocStatKS))) + logVolUnitBall)
#if MAHAL_TRANSFORM_ENABLED
            beta = exp(log(beta) - sum(log(choDia)))
#endif
            alpha = real(minLocStatKS,RK)
            densityAvg = alpha / beta
            densityStd = sqrt(alpha) / beta
            volumeAvg = np * beta / (alpha - 1._RK) ! effSamSize / densityAvg
            volumeStd = np * beta / (alpha - 1._RK) / sqrt(alpha - 2._RK) ! effSamSize / densityStd
            write(*, "(*(g0.8,:,', '))")
            write(*, "(*(g0.8,:,', '))") "MinAdjustedStatKS"
            write(*, "(*(g0.8,:,', '))") "alpha, beta", alpha, beta
            write(*, "(*(g0.8,:,', '))") "densityAvg, densityStd", densityAvg, densityStd
            write(*, "(*(g0.8,:,', '))") "volumeAvg, volumeStd", volumeAvg, volumeStd
            write(*, "(*(g0.8,:,', '))") "(volumeAvg - volumeRef) / volumeStd", (volumeAvg - volumeRef) / volumeStd
        end block

        read(*,*)
#endif

        !test%File = test%openFile() ! label = "LogVolEstimate")
        !write(test%File%unit, "(*(g0,:,','))") ClusteredPoint%LogDensity, logDensityFromExpectedMinPairDist
        !close(test%File%unit)

    end function test_getPosteriorDistSortedExpDiff_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_processPoisson ! LCOV_EXCL_LINE