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

!>  \brief This module contains tests of the module [pm_clustering](@ref pm_clustering).
!>  \author Amir Shahmoradi

module test_pm_clustering

    use pm_test, only: test_type, LK
    use pm_kind, only: LK
    use pm_clustering
    implicit none

    private
    public :: setTest
    type(test_type) :: test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 1000
        real(RK), allocatable   :: Point(:,:)
    contains
        procedure, pass :: read => readTestData
    end type TestData_type

    type(TestData_type) :: TestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call TestData%read()
        call test%run(test_benchmark_1, SK_"test_benchmark_1")
        call test%run(test_runKmeans_1, SK_"test_runKmeans_1")
        call test%run(test_runKmeans_2, SK_"test_runKmeans_2")
        call test%run(test_runKmeans_3, SK_"test_runKmeans_3")
        call test%run(test_runKmeans_4, SK_"test_runKmeans_4")
        call test%run(test_setKmeans_1, SK_"test_setKmeans_1")
        call test%run(test_setKmeans_2, SK_"test_setKmeans_2")
        call test%run(test_setKmeans_3, SK_"test_setKmeans_3")
        call test%run(test_setKmeans_4, SK_"test_setKmeans_4")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine readTestData(TestData)
        implicit none
        class(TestData_type), intent(inout) :: TestData
        integer(IK) :: fileUnit, ip
        open( file = test%dir%inp//"/test_pm_clustering@points.txt" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
#if         INTEL_ENABLED && WINDOWS_ENABLED
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        if (allocated(TestData%Point)) deallocate(TestData%Point) ! LCOV_EXCL_LINE
        allocate(TestData%Point(TestData%nd,TestData%np))
        read(fileUnit,*)
        do ip = 1, TestData%np
            read(fileUnit,*) TestData%Point(1:TestData%nd,ip)
        end do
        close(fileUnit)
    end subroutine readTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_runKmeans_1() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK), parameter      :: nc = 3_IK
        logical(LK)                 :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true._LK

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            )

        ! write results to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. .not. Kmeans%err%occurred
        assertion = assertion .and. Kmeans%err%stat /= 1_IK
        assertion = assertion .and. Kmeans%err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `setKmeans()` by passing a fixed initial set of cluster centers to the Kmeans constructor.
    function test_runKmeans_2() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        real(RK)    , allocatable   :: InitCenter(:,:)
        logical(LK)                 :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true._LK

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. .not. Kmeans%err%occurred
        assertion = assertion .and. Kmeans%err%stat /= 1_IK
        assertion = assertion .and. Kmeans%err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> If the optional input argument `niterMax` is specified, the output value for `niter` must not go beyond in the input value.
    !> In addition, if the specified value for `niterMax` has reached, the procedure must return with error stat code of `1`.
    function test_runKmeans_3() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: niterMax = 0_IK
        logical(LK)                 :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true._LK

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , niterMax = niterMax & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. Kmeans%err%occurred
        assertion = assertion .and. Kmeans%err%stat == 1_IK
        assertion = assertion .and. Kmeans%niter == niterMax
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The function `setKmeans()` must function properly for reasonable optional input values of `nfailMax` and `relTol`.
    function test_runKmeans_4() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nfailMax = 100_IK
        real(RK)    , parameter     :: relTol = 1.e-8_RK
        logical(LK)                 :: assertion
        type(Kmeans_type)           :: Kmeans

        assertion = .true._LK

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                            , relTol = relTol & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. .not. Kmeans%err%occurred
        assertion = assertion .and. Kmeans%err%stat /= 1_IK
        assertion = assertion .and. Kmeans%err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `setKmeans()` by passing a number of tries to find the more optimal Kmeans clustering.
    function test_setKmeans_1() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nt = 2_IK + nint(log(real(nc)))
        real(RK)    , allocatable   :: InitCenter(:,:)
        logical(LK)                 :: assertion
        type(Kmeans_type)           :: Kmeans

        assertion = .true._LK

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. .not. Kmeans%err%occurred
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. Kmeans%err%stat /= 1_IK .and. Kmeans%err%stat /= 2_IK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_setKmeans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The component `index` must be properly set by [pm_clustering::setKmeans](@ref pm_clustering::setKmeans) when it is given as input.
    function test_setKmeans_2() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nt = 2_IK + nint(log(real(nc)))
        Integer(IK) , allocatable   :: PointIndex(:)
        real(RK)    , allocatable   :: InitCenter(:,:), Point(:,:)
        logical(LK)                 :: assertion
        type(Kmeans_type)           :: Kmeans
        integer(IK)                 :: ip, ic

        assertion = .true._LK

        Point = TestData%Point
        PointIndex = [(ip,ip=1,TestData%np)]
        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )

        assertion = assertion .and. .not. Kmeans%err%occurred
        if (.not. assertion) return

        call Kmeans%getProp(nd = TestData%nd, np = TestData%np, Point = Point, index = PointIndex)

        assertion = assertion .and. .not. Kmeans%err%occurred
        call test%assert(assertion)

        ! write data to output for further investigation

        test%File = test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        close(test%File%unit)

        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. Kmeans%err%stat /= 1_IK .and. Kmeans%err%stat /= 2_IK
        assertion = assertion .and. all(Kmeans%membership > 0_IK) .and. all(Kmeans%membership < nc + 1)
        assertion = assertion .and. all(Kmeans%minDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np
        assertion = assertion .and. all(TestData%Point(:,PointIndex) == Point)
        do ic = 1, nc
            assertion = assertion .and. all( Kmeans%membership(Kmeans%Prop%cumSumSize(ic-1)+1:Kmeans%Prop%cumSumSize(ic)) == Kmeans%membership(Kmeans%Prop%cumSumSize(ic)) )
        end do
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0.15,:,' '))")
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership, mask = Kmeans%membership < 1_IK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership, mask = Kmeans%membership > nc)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq, mask = Kmeans%minDistanceSq < 0._RK)
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
            write(test%disp%unit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_setKmeans_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The component `index` must be properly set by [pm_clustering::setKmeans](@ref pm_clustering::setKmeans) when it is given as input.
    function test_setKmeans_3() result(assertion)
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: nc = 2_IK
        integer(IK) , parameter     :: nd = 2_IK
        integer(IK) , parameter     :: np = 11_IK
        integer(IK) , parameter     :: nt = 10_IK + 2_IK + nint(log(real(nc)))
        real(RK)    , parameter     :: LOG_PI = log(acos(-1._RK))
        real(RK)    , parameter     :: pointLogVolNormed = log(1.e-1_RK) - LOG_PI
        integer(IK) , parameter     :: minSize = 1_IK
        Integer(IK)                 :: PointIndex(np)
        real(RK)                    :: Point_ref(nd,np)
        real(RK)                    :: Point(nd,np)
        logical(LK)                 :: assertion
        type(Kmeans_type)           :: Kmeans
        integer(IK)                 :: ip

        assertion = .true._LK

        Point_ref = reshape([ 0.126986816293506_RK, 0.957166948242946_RK &
                            , 0.913375856139019_RK, 0.485375648722841_RK &
                            , 0.632359246225410_RK, 0.800280468888800_RK &
                            , 0.097540404999410_RK, 0.141886338627215_RK &
                            , 0.278498218867048_RK, 0.421761282626275_RK &
                            , 0.546881519204984_RK, 0.915735525189067_RK &
                            , 0.957506835434298_RK, 0.792207329559554_RK &
                            , 0.964888535199277_RK, 0.959492426392903_RK &
                            , 0.157613081677548_RK, 0.655740699156587_RK &
                            , 0.970592781760616_RK, 0.035711678574190_RK &
                            , 10.00000000000000_RK, 10.00000000000000_RK ], shape = shape(Point_ref))

        Point = Point_ref
        PointIndex = [(ip,ip=1,np)]

        Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , minSize = minSize & ! LCOV_EXCL_LINE
                            )

        assertion = assertion .and. (.not. Kmeans%err%occurred .or. (Kmeans%err%occurred .and. Kmeans%err%stat == 2_IK))
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "No error must occur, otherwise stat must be 2 indicating zero-sized cluster."
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all(Kmeans%membership(1:10) == Kmeans%membership(1))
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "The first 10 memberships must be equal."
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership    =", Kmeans%membership
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        call Kmeans%getProp(nd = nd, np = np, Point = Point, index = PointIndex, pointLogVolNormed = pointLogVolNormed)

        assertion = assertion .and. all(Point_ref(:,PointIndex) == Point)
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point_ref(:,PointIndex) == Point"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point_ref(:,PointIndex)  =", Point_ref(:,PointIndex)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point                    =", Point
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all(Kmeans%Prop%EffectiveSize == Kmeans%Size)
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "all(Kmeans%Prop%EffectiveSize == Kmeans%Size)"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%EffectiveSize    =", Kmeans%Prop%EffectiveSize
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size                  =", Kmeans%Size
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. ( abs(Kmeans%Prop%LogVolNormed(2) - pointLogVolNormed) < 1.e-10_RK .or. abs(Kmeans%Prop%LogVolNormed(1) - pointLogVolNormed) < 1.e-10_RK )
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogVolNormed(2)  == pointLogVolNormed"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogVolNormed(2)  =", Kmeans%Prop%LogVolNormed(2)
                write(test%disp%unit,"(*(g0.15,:,' '))") "pointLogVolNormed            =", pointLogVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        !assertion = assertion .and. ( abs(Kmeans%Prop%LogDenNormed(2) + pointLogVolNormed) < 1.e-10_RK .or. abs(Kmeans%Prop%LogDenNormed(1) + pointLogVolNormed) < 1.e-10_RK )
        !if (.not. assertion) then
        !    ! LCOV_EXCL_START
        !    if (test%traceable) then
        !        write(test%disp%unit,"(*(g0.15,:,' '))")
        !        write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogDenNormed(2)  == pointLogVolNormed"
        !        write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogDenNormed(2)  =", Kmeans%Prop%LogDenNormed(2)
        !        write(test%disp%unit,"(*(g0.15,:,' '))") "-pointLogVolNormed           =", -pointLogVolNormed
        !        write(test%disp%unit,"(*(g0.15,:,' '))")
        !    end if
        !    ! LCOV_EXCL_STOP
        !    return
        !end if

        assertion = assertion .and. abs(Kmeans%Prop%logSumVolNormed + 0.552093409710310_RK) < 1.e-10_RK
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%logSumVolNormed  == -0.552093409710310_RK"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%logSumVolNormed  =", Kmeans%Prop%logSumVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

    end function test_setKmeans_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> When the `pointLogVolNormed` is missing, the properties of singular clusters must be correctly computed from 
    !> the properties of non-singular clusters.
    function test_setKmeans_4() result(assertion)
        use pm_mathCumSum, only: getCumSum
        use pm_val2str, only: getStr
        use pm_kind, only: IK, RK
        implicit none
        integer(IK) , parameter     :: nc = 2_IK
        integer(IK) , parameter     :: nd = 2_IK
        integer(IK) , parameter     :: np = 11_IK
        integer(IK) , parameter     :: nt = 10_IK + 2_IK + nint(log(real(nc)))
        real(RK)    , parameter     :: logSumVolNormed = -0.513657091888111_RK
        real(RK)    , parameter     :: LOG_PI = log(acos(-1._RK))
        real(RK)    , parameter     :: pointLogVolNormed = log(1.e-1_RK) - LOG_PI
        real(RK)    , parameter     :: LogVolNormed(nc) = [-0.608967271692436_RK, -2.91155236468648_RK]
       !real(RK)    , parameter     :: LogDenNormed(nc) = [2.91155236468648_RK, 2.91155236468648_RK]
        real(RK)    , parameter     :: choLowCovUpp(nd,nd,nc) = reshape([ 0.140347692021171_RK &
                                                                        , 0.273777628663980E-1_RK &
                                                                        , 0.492129448979672E-2_RK &
                                                                        , 0.111901976129643_RK &
                                                                        , 0.543912292747094E-1_RK &
                                                                        , 0.00000000000000_RK &
                                                                        , 0.00000000000000_RK &
                                                                        , 0.543912292747094E-1_RK &
                                                                        ], shape = shape(choLowCovUpp))
        real(RK)    , parameter     :: invCov(nd,nd,nc) = reshape(   [ 1.64294297527196_RK &
                                                                        , -0.722543648549051E-1_RK &
                                                                        , -0.722543648549051E-1_RK &
                                                                        , 2.06058250870097_RK &
                                                                        , 18.3853171427581_RK &
                                                                        , 0.00000000000000_RK &
                                                                        , 0.00000000000000_RK &
                                                                        , 18.3853171427581_RK &
                                                                        ], shape = shape(invCov))
        real(RK)    , parameter     :: choDia(nd,nc) = reshape( [ 0.780771367974080_RK &
                                                                , 0.696634527157956_RK &
                                                                , 0.233219272948677_RK &
                                                                , 0.233219272948677_RK &
                                                                ], shape = shape(choDia))
        integer(IK) , parameter     :: minSize = 1_IK
        Integer(IK)                 :: PointIndex(np)
        real(RK)                    :: Point_ref(nd,np)
        real(RK)                    :: Point(nd,np)
        logical(LK)                 :: assertion
        type(Kmeans_type)           :: Kmeans
        integer(IK)                 :: ip

        assertion = .true._LK

        Point_ref = reshape([ 0.126986816293506_RK, 0.957166948242946_RK &
                            , 0.913375856139019_RK, 0.485375648722841_RK &
                            , 0.632359246225410_RK, 0.800280468888800_RK &
                            , 0.097540404999410_RK, 0.141886338627215_RK &
                            , 0.278498218867048_RK, 0.421761282626275_RK &
                            , 0.546881519204984_RK, 0.915735525189067_RK &
                            , 0.957506835434298_RK, 0.792207329559554_RK &
                            , 0.964888535199277_RK, 0.959492426392903_RK &
                            , 0.157613081677548_RK, 0.655740699156587_RK &
                            , 0.970592781760616_RK, 0.035711678574190_RK &
                            , 10.00000000000000_RK, 10.00000000000000_RK ], shape = shape(Point_ref))

        Point = Point_ref
        PointIndex = [(ip,ip=1,np)]

        Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , minSize = minSize & ! LCOV_EXCL_LINE
                            )

        assertion = assertion .and. (.not. Kmeans%err%occurred .or. (Kmeans%err%occurred .and. Kmeans%err%stat == 2_IK))
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "No error must occur, otherwise stat must be 2 indicating zero-sized cluster."
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        ! Test getProp() with no input pointLogVolNormed or index

        call Kmeans%getProp(nd = nd, np = np, Point = Point, inclusionFraction = 0._RK)

        assertion = assertion .and. all(Point_ref(:,Kmeans%Prop%index) == Point)
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point_ref(:,Kmeans%Prop%index) == Point"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point_ref    =", Point_ref(:,Kmeans%Prop%index)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Point        =", Point
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. abs(Kmeans%Prop%logSumVolNormed - logSumVolNormed) < 1.e-8_RK
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "all(Kmeans%Prop%logSumVolNormed == logSumVolNormed"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%logSumVolNormed      =", Kmeans%Prop%logSumVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))") "logSumVolNormed                  =", logSumVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. all(Kmeans%Prop%EffectiveSize == Kmeans%Size)
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "all(Kmeans%Prop%EffectiveSize == Kmeans%Size)"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%EffectiveSize    =", Kmeans%Prop%EffectiveSize
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%Size             =", Kmeans%Size
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. ( all(abs(Kmeans%Prop%LogVolNormed - LogVolNormed) < 1.e-10_RK) .or. all(abs(Kmeans%Prop%LogVolNormed(nc:1:-1) - LogVolNormed) < 1.e-10_RK) )
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogVolNormed(nc:1:-1)    == LogVolNormed)"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogVolNormed             =", Kmeans%Prop%LogVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))") "LogVolNormed                         =", LogVolNormed
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        !assertion = assertion .and. ( all(abs(Kmeans%Prop%LogDenNormed - LogDenNormed) < 1.e-10_RK) .or. all(abs(Kmeans%Prop%LogDenNormed(nc:1:-1) - LogDenNormed) < 1.e-10_RK) )
        !if (.not. assertion) then
        !    ! LCOV_EXCL_START
        !    if (test%traceable) then
        !        write(test%disp%unit,"(*(g0.15,:,' '))")
        !        write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogDenNormed(nc:1:-1)    == LogDenNormed)"
        !        write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%LogDenNormed             =", Kmeans%Prop%LogDenNormed
        !        write(test%disp%unit,"(*(g0.15,:,' '))")
        !    end if
        !    ! LCOV_EXCL_STOP
        !    return
        !end if

        assertion = assertion .and. all(Kmeans%Prop%cumSumSize(1:) == getCumSum(Kmeans%Size(1:nc))) .and. Kmeans%Prop%cumSumSize(0) == 0_IK
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "all(Kmeans%Prop%cumSumSize == [0_IK,10_IK,11_IK])"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%cumSumSize     =", Kmeans%Prop%cumSumSize
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. ( all(abs(Kmeans%Prop%choLowCovUpp - choLowCovUpp) < 1.e-10_RK) .or. all(abs(Kmeans%Prop%choLowCovUpp(:,:,nc:1:-1) - choLowCovUpp) < 1.e-10_RK) )
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%choLowCovUpp(:,:,nc:1:-1)    == choLowCovUpp"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%choLowCovUpp                 =", Kmeans%Prop%choLowCovUpp
                write(test%disp%unit,"(*(g0.15,:,' '))") "choLowCovUpp                             =", choLowCovUpp
                write(test%disp%unit,"(*(g0.15,:,' '))") "difference                               =", abs(Kmeans%Prop%choLowCovUpp - choLowCovUpp)
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. ( all(abs(Kmeans%Prop%invCov - invCov) < 1.e-10_RK) .or. all(abs(Kmeans%Prop%invCov(:,:,nc:1:-1) - invCov) < 1.e-10_RK) )
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%invCov(:,:,nc:1:-1)   == invCov"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%invCov                =", Kmeans%Prop%invCov
                write(test%disp%unit,"(*(g0.15,:,' '))") "invCov                            =", invCov
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        assertion = assertion .and. ( all(abs(Kmeans%Prop%choDia - choDia) < 1.e-10_RK) .or. all(abs(Kmeans%Prop%choDia(:,nc:1:-1) - choDia) < 1.e-10_RK) )
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (test%traceable) then
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%choDia(:,nc:1:-1)    == choDia"
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%choDia               =", Kmeans%Prop%choDia
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Prop%choDia               =", choDia
                write(test%disp%unit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end if
        call test%assert(assertion)

        ! write data to output for further investigation

        !test%File = test%openFile()
        !call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = test%File%unit)
        !close(test%File%unit)

    end function test_setKmeans_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Calling the Kmeans routine repeatedly should not cause any errors.
    !>  This test is also used for benchmarking the performances of different implementations of the Kmeans algorithm.
    function test_benchmark_1() result(assertion)
        use pm_distUnif, only: getUnifRand
        use pm_kind, only: IK, RK
        use pm_val2str, only: getStr
        implicit none
        integer(IK) , parameter     :: ncMax = 3_IK
        integer(IK) , parameter     :: nfailMax = 100_IK
        real(RK)    , parameter     :: relTol = 1.e-8_RK
        logical(LK)                 :: assertion

        type(Kmeans_type)           :: Kmeans
        integer(IK)                 :: i, np, nc

        assertion = .true._LK

        do i = 1_IK, 100_IK

            nc = ncMax ! getUnifRand(2, ncMax)
            np = TestData%np ! getUnifRand(100,TestData%np)
            Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = 1_IK & ! LCOV_EXCL_LINE
                                , Point = TestData%Point(1:TestData%nd,1:np)& ! LCOV_EXCL_LINE
                                , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                                , relTol = relTol & ! LCOV_EXCL_LINE
                                )

            assertion = assertion .and. .not. Kmeans%err%occurred
            assertion = assertion .and. Kmeans%err%stat /= 1_IK
            assertion = assertion .and. Kmeans%err%stat /= 2_IK
            assertion = assertion .and. Kmeans%potential > 0._RK
            assertion = assertion .and. all(Kmeans%membership(1:np) > 0_IK) .and. all(Kmeans%membership(1:np) < nc + 1)
            assertion = assertion .and. all(Kmeans%minDistanceSq(1:np) > 0_IK)
            assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

            if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0.15,:,' '))")
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership < 1    =", pack(Kmeans%membership(1:np), mask = Kmeans%membership(1:np) < 1_IK)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%membership > nc   =", pack(Kmeans%membership(1:np), mask = Kmeans%membership(1:np) > nc)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%minDistanceSq < 0 =", pack(Kmeans%minDistanceSq(1:np), mask = Kmeans%minDistanceSq(1:np) < 0._RK)
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%occurred      =", Kmeans%err%occurred
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%err%stat          =", Kmeans%err%stat
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
                write(test%disp%unit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
                write(test%disp%unit,"(*(g0.15,:,' '))")
                return
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_benchmark_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_clustering ! LCOV_EXCL_LINE