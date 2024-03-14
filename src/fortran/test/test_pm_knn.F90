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

!>  \brief This module contains tests of the module [pm_knn](@ref pm_knn).
!>  \author Amir Shahmoradi

module test_pm_knn

    use pm_test, only: test_type, LK

    use pm_container, only: IV => cvi_pdt, RV => cvr_pdt
    use pm_kind, only: LK
    use pm_knn

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

    type(TestData_type) :: TestData

    interface TestData_type
        module procedure :: contructTestData
    end interface TestData_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        TestData = TestData_type()
        call test%run(test_getHub_1, SK_"test_getHub_1")
        call test%run(test_getEucDistSq_1, SK_"test_getEucDistSq_1")
        call test%run(test_getPairDistSq_1, SK_"test_getPairDistSq_1")
        call test%run(test_getDistSortedExpDiff_1, SK_"test_getDistSortedExpDiff_1")
        call test%run(test_getDistSortedExpDiff_2, SK_"test_getDistSortedExpDiff_2")
        call test%run(test_getDistSortedExpDiff_3, SK_"test_getDistSortedExpDiff_3")
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

    function test_getEucDistSq_1() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)         :: assertion
        real(RK), parameter :: Point1(*) = [0._RK, 1._RK, 2._RK, 3._RK, 4._RK]
        real(RK), parameter :: Point2(*) = Point1 + 1._RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: distanceSq_ref = norm2(Point2-Point1)**2
        real(RK)            :: distanceSq
        real(RK)            :: difference
        distanceSq = getDisEuclid(point1, point2, euclidsq)
        difference = abs(distanceSq - distanceSq_ref) / distanceSq_ref
        assertion  = difference < tolerance
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "distanceSq_ref    = ", distanceSq_ref
            write(test%disp%unit,"(*(g0,:,' '))") "distanceSq        = ", distanceSq
            write(test%disp%unit,"(*(g0,:,' '))") "difference        = ", difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEucDistSq_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getPairDistSq()`.
    function test_getPairDistSq_1() result(assertion)
        use pm_distanceEuclid, only: getDisMatEuclid, euclidsq
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)                 :: assertion
        real(RK), allocatable       :: Dist(:,:)
        real(RK), allocatable       :: diff(:,:)
        real(RK), parameter         :: tolerance = 1.e-10_RK
        integer(IK)                 :: ip, jp

        assertion = .true._LK

        Dist = getDisMatEuclid(TestData%Point, diag = 0._RK)
        diff = abs( (Dist - TestData%Dist) )

        assertion = all( diff < tolerance )

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Distance, Distance_ref, diff"
            do jp = 1, TestData%np
                do ip = 1, TestData%np
                    write(test%disp%unit,"(*(g0.15,:,' '))") Dist(ip,jp), TestData%Dist(ip,jp), diff(ip,jp)
                end do
            end do
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPairDistSq_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `setDisSortedExpDiff()` for an even value of `nd`.
    function test_getDistSortedExpDiff_1() result(assertion)
        use pm_kind, only: IK, RK
        use pm_err, only: err_type
        implicit none
        logical(LK)                 :: assertion
        real(RK), allocatable       :: diff(:)
        real(RK), allocatable       :: RefPoint(:)
        real(RK), allocatable       :: disSortedExpDiff(:)
        real(RK), parameter         :: DistSortedExpDiff_ref(*) =   [ 0.05518433594135580_RK & ! LCOV_EXCL_LINE
                                                                    , 0.06693027963675390_RK & ! LCOV_EXCL_LINE
                                                                    , 0.01936935163401820_RK & ! LCOV_EXCL_LINE
                                                                    , 0.03354993653506870_RK & ! LCOV_EXCL_LINE
                                                                    , 0.01353010700222810_RK & ! LCOV_EXCL_LINE
                                                                    , 0.00611135862293311_RK & ! LCOV_EXCL_LINE
                                                                    , 0.10002225854507400_RK & ! LCOV_EXCL_LINE
                                                                    , 0.03619526286475090_RK & ! LCOV_EXCL_LINE
                                                                    , 0.09636174928998410_RK & ! LCOV_EXCL_LINE
                                                                    , 0.00976657158542438_RK & ! LCOV_EXCL_LINE
                                                                    ]
        real(RK), parameter         :: tolerance = 1.e-10_RK
        integer(IK), allocatable    :: index(:)
        type(err_type)              :: Err
        integer(IK)                 :: ip

        assertion = .true._LK

        if (allocated(RefPoint)) deallocate(RefPoint) ! LCOV_EXCL_LINE
        allocate(RefPoint(TestData%nd), source = 0.5_RK)
        if (allocated(disSortedExpDiff)) deallocate(disSortedExpDiff) ! LCOV_EXCL_LINE
        allocate(disSortedExpDiff(TestData%np))
        if (allocated(index)) deallocate(index) ! LCOV_EXCL_LINE
        allocate(index(TestData%np))

        call setDisSortedExpDiff   ( nd = TestData%nd & ! LCOV_EXCL_LINE
                                    , np = TestData%np & ! LCOV_EXCL_LINE
                                    , Point = TestData%Point & ! LCOV_EXCL_LINE
                                    , RefPoint = RefPoint & ! LCOV_EXCL_LINE
                                    , disSortedExpDiff = disSortedExpDiff & ! LCOV_EXCL_LINE
                                    , index = index & ! LCOV_EXCL_LINE
                                    )
        assertion = .not. err%occurred
        call test%assert(assertion)

        diff = abs( (disSortedExpDiff - DistSortedExpDiff_ref) )
        assertion = all( diff < tolerance )
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "disSortedExpDiff, DistSortedExpDiff_ref, diff"
            do ip = 1, TestData%np
                write(test%disp%unit,"(*(g0.15,:,' '))") disSortedExpDiff(ip), DistSortedExpDiff_ref(ip), diff(ip)
            end do
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDistSortedExpDiff_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `setDisSortedExpDiff()` for an even value of `nd` but with a reference point from within input set of points.
    function test_getDistSortedExpDiff_2() result(assertion)
        use pm_kind, only: IK, RK
        use pm_err, only: err_type
        implicit none
        logical(LK)                 :: assertion
        real(RK), allocatable       :: diff(:)
        real(RK), allocatable       :: RefPoint(:)
        real(RK), allocatable       :: disSortedExpDiff(:)
        real(RK), parameter         :: DistSortedExpDiff_ref(*) =   [ 0.0000000000000000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0693595837531895_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0621850083406652_RK & ! LCOV_EXCL_LINE
                                                                    , 0.1736352649921330_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0108602911297795_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0222521392316455_RK & ! LCOV_EXCL_LINE
                                                                    , 0.2599906876136520_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0297462057145467_RK & ! LCOV_EXCL_LINE
                                                                    , 0.0152054907472889_RK & ! LCOV_EXCL_LINE
                                                                    , 0.1170517778083240_RK & ! LCOV_EXCL_LINE
                                                                    ]
        real(RK), parameter         :: tolerance = 1.e-10_RK
        integer(IK), allocatable    :: index(:)
        type(err_type)              :: Err
        integer(IK)                 :: ip

        assertion = .true._LK

        if (allocated(RefPoint)) deallocate(RefPoint) ! LCOV_EXCL_LINE
        allocate(RefPoint(TestData%nd), source = TestData%Point(1:TestData%nd,1))
        if (allocated(disSortedExpDiff)) deallocate(disSortedExpDiff) ! LCOV_EXCL_LINE
        allocate(disSortedExpDiff(TestData%np))
        if (allocated(index)) deallocate(index) ! LCOV_EXCL_LINE
        allocate(index(TestData%np))

        call setDisSortedExpDiff   ( nd = TestData%nd & ! LCOV_EXCL_LINE
                                    , np = TestData%np & ! LCOV_EXCL_LINE
                                    , Point = TestData%Point & ! LCOV_EXCL_LINE
                                    , RefPoint = RefPoint & ! LCOV_EXCL_LINE
                                    , disSortedExpDiff = disSortedExpDiff & ! LCOV_EXCL_LINE
                                    , index = index & ! LCOV_EXCL_LINE
                                    )
        assertion = .not. err%occurred
        call test%assert(assertion)

        diff = abs( (disSortedExpDiff - DistSortedExpDiff_ref) )
        assertion = all( diff < tolerance )
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "disSortedExpDiff, DistSortedExpDiff_ref, diff"
            do ip = 1, TestData%np
                write(test%disp%unit,"(*(g0.15,:,' '))") disSortedExpDiff(ip), DistSortedExpDiff_ref(ip), diff(ip)
            end do
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDistSortedExpDiff_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `setDisSortedExpDiff()` for an odd value of `nd` with a reference point from within input set of points.
    function test_getDistSortedExpDiff_3() result(assertion)
        use pm_kind, only: IK, RK
        use pm_err, only: err_type
        implicit none
        logical(LK)                 :: assertion
        real(RK), allocatable       :: diff(:)
        real(RK), allocatable       :: Point(:)
        real(RK), allocatable       :: RefPoint(:)
        real(RK), allocatable       :: disSortedExpDiff(:)
        real(RK), parameter         :: DistSortedExpDiff_ref(*) =   [ 0.000000000000000000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.120885137189500000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.015726743050333000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.070265549615960000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.061505870482143000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.253398949683816000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.156886479354146000_RK & ! LCOV_EXCL_LINE
                                                                    , 0.000339887191352029_RK & ! LCOV_EXCL_LINE
                                                                    , 0.007381699764978930_RK & ! LCOV_EXCL_LINE
                                                                    , 0.005704246561339060_RK & ! LCOV_EXCL_LINE
                                                                    ]
        real(RK), parameter         :: tolerance = 1.e-10_RK
        integer(IK), allocatable    :: index(:)
        type(err_type)              :: Err
        integer(IK)                 :: ip

        assertion = .true._LK

        if (allocated(RefPoint)) deallocate(RefPoint) ! LCOV_EXCL_LINE
        allocate(RefPoint(TestData%nd), source = TestData%Point(1,1))
        if (allocated(disSortedExpDiff)) deallocate(disSortedExpDiff) ! LCOV_EXCL_LINE
        allocate(disSortedExpDiff(TestData%np))
        if (allocated(index)) deallocate(index) ! LCOV_EXCL_LINE
        allocate(index(TestData%np))

        Point = TestData%Point(1,:)
        call setDisSortedExpDiff   ( nd = 1_IK & ! LCOV_EXCL_LINE
                                    , np = TestData%np & ! LCOV_EXCL_LINE
                                    , Point = Point & ! LCOV_EXCL_LINE
                                    , RefPoint = RefPoint & ! LCOV_EXCL_LINE
                                    , disSortedExpDiff = disSortedExpDiff & ! LCOV_EXCL_LINE
                                    , index = index & ! LCOV_EXCL_LINE
                                    )
        assertion = .not. err%occurred
        call test%assert(assertion)

        diff = abs( (disSortedExpDiff - DistSortedExpDiff_ref) )
        assertion = all( diff < tolerance )
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "disSortedExpDiff, DistSortedExpDiff_ref, diff"
            do ip = 1, TestData%np
                write(test%disp%unit,"(*(g0.15,:,' '))") disSortedExpDiff(ip), DistSortedExpDiff_ref(ip), diff(ip)
            end do
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDistSortedExpDiff_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getHub()`.
    function test_getHub_1() result(assertion)
        use pm_distanceEuclid, only: getDisMatEuclid, euclidsq
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)                 :: assertion
        real(RK), parameter         :: tolerance = 1.e-10_RK
        type(hub_type)              :: Hub
        integer(IK)                 :: i

        assertion = .true._LK

        Hub = hub_type  ( np = TestData%np & ! LCOV_EXCL_LINE
                        , pairDisSq = getDisMatEuclid(TestData%Point, diag = 0._RK, method = euclidsq) & ! LCOV_EXCL_LINE
                        )

        assertion = .not. Hub%err%occurred
        call test%assert(assertion)

        assertion = assertion .and. all( Hub%minDisEdge%index%val == TestData%HubMinDistEdgeIndex )
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%minDisEdge%index        ", Hub%minDisEdge%index%val
            write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubMinDistEdgeIndex ", TestData%HubMinDistEdgeIndex
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all( abs((Hub%minDisEdge%LenSq%val - TestData%HubMinDistEdgeLenSq)/TestData%HubMinDistEdgeLenSq) < tolerance )
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%minDisEdge%LenSq        ", Hub%minDisEdge%LenSq%val
            write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubMinDistEdgeLenSq ", TestData%HubMinDistEdgeLenSq
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all( Hub%NodeIndex == TestData%HubNodeIndex )
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%NodeIndex        ", Hub%NodeIndex
            write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubNodeIndex", TestData%HubNodeIndex
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all( Hub%EdgeCount == TestData%HubEdgeCount )
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%EdgeCount        ", Hub%EdgeCount
            write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubEdgeCount", TestData%HubEdgeCount
            write(test%disp%unit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        do i = 1, Hub%nh
            assertion = assertion .and. all(Hub%EdgeIndex(i)%val == TestData%HubEdgeIndex(i)%val)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%EdgeIndex(i)%val           ", Hub%EdgeIndex(i)%val
                write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubEdgeIndex(i)%val   ", TestData%HubEdgeIndex(i)%val
                write(test%disp%unit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
            if (.not. assertion) exit ! LCOV_EXCL_LINE
        end do
        call test%assert(assertion)

        do i = 1, Hub%nh
            assertion = assertion .and. all( abs((Hub%EdgeLenSq(i)%val - TestData%HubEdgeLenSq(i)%val)/TestData%HubEdgeLenSq(i)%val) < tolerance )
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Hub%EdgeLenSq(i)%val           ", Hub%EdgeLenSq(i)%val
                write(test%disp%unit,"(*(g0.15,:,' '))") "TestData%HubEdgeLenSq(i)%val   ", TestData%HubEdgeLenSq(i)%val
                write(test%disp%unit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
            if (.not. assertion) exit ! LCOV_EXCL_LINE
        end do

    end function test_getHub_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_knn ! LCOV_EXCL_LINE