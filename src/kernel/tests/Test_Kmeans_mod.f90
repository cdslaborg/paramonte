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

!>  \brief This module contains tests of the module [Kmeans_mod](@ref Kmeans_mod).
!>  \author Amir Shahmoradi

module Test_Kmeans_mod

    use Test_mod, only: Test_type
    use Kmeans_mod
    implicit none

    private
    public :: test_Kmeans

    type(Test_type) :: Test

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

    subroutine test_Kmeans()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call TestData%read()
        call Test%run(test_benchmark_1, "test_benchmark_1")
        call Test%run(test_runKmeans_1, "test_runKmeans_1")
        call Test%run(test_runKmeans_2, "test_runKmeans_2")
        call Test%run(test_runKmeans_3, "test_runKmeans_3")
        call Test%run(test_runKmeans_4, "test_runKmeans_4")
        call Test%run(test_getKmeans_1, "test_getKmeans_1")
        call Test%run(test_getKmeans_2, "test_getKmeans_2")
        call Test%finalize()
    end subroutine test_Kmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine readTestData(TestData)
        implicit none
        class(TestData_type), intent(inout) :: TestData
        integer(IK) :: fileUnit, ip
        open( file = Test%inDir//"/Test_Kmeans_mod@points.txt" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
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
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK), parameter      :: nc = 3_IK
        logical                     :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true.

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            )

        ! write results to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. .not. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runKmeans()` by passing a fixed initial set of cluster centers to the Kmeans constructor.
    function test_runKmeans_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        real(RK)    , allocatable   :: InitCenter(:,:)
        logical                     :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true.

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. .not. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> If the optional input argument `niterMax` is specified, the output value for `niter` must not go beyond in the input value.
    !> In addition, if the specified value for `niterMax` has reached, the procedure must return with error stat code of `1`.
    function test_runKmeans_3() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: niterMax = 1_IK
        logical                     :: assertion

        type(Kmeans_type)           :: Kmeans

        assertion = .true.

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , niterMax = niterMax & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%Err%stat == 1_IK
        assertion = assertion .and. Kmeans%niter == niterMax
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The function `runKmeans()` must function properly for reasonable optional input values of `nfailMax` and `relTol`.
    function test_runKmeans_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nfailMax = 100_IK
        real(RK)    , parameter     :: relTol = 1.e-8_RK
        logical                     :: assertion
        type(Kmeans_type)           :: Kmeans

        assertion = .true.

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = 1_IK & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                            , relTol = relTol & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. .not. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `getKmeans()` by passing a number of tries to find the more optimal Kmeans clustering.
    function test_getKmeans_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nt = 2_IK + nint(log(real(nc)))
        real(RK)    , allocatable   :: InitCenter(:,:)
        logical                     :: assertion
        type(Kmeans_type)           :: Kmeans

        assertion = .true.

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            , propEnabled = .true. & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. .not. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. Kmeans%Err%stat /= 1_IK .and. Kmeans%Err%stat /= 2_IK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getKmeans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The component `Index` must be properly set by [getKmeans](@ref kmeans_mod::gekmeans) when it is given as input.
    function test_getKmeans_2() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nt = 2_IK + nint(log(real(nc)))
        Integer(IK) , allocatable   :: Index(:)
        real(RK)    , allocatable   :: InitCenter(:,:), Point(:,:)
        logical                     :: assertion
        type(Kmeans_type)           :: Kmeans
        integer(IK)                 :: ip, ic

        assertion = .true.

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])
        Point = TestData%Point
        Index = [(ip,ip=1,TestData%np)]

        Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Index = Index & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            , propEnabled = .true. & ! LCOV_EXCL_LINE
                            )

        ! write data to output for further investigation

        Test%File = Test%openFile()
        call writeKmeans(Kmeans = Kmeans, Point = TestData%Point, nd = TestData%nd, np = TestData%np, fileUnit = Test%File%unit)
        close(Test%File%unit)

        assertion = assertion .and. .not. Kmeans%Err%occurred
        assertion = assertion .and. Kmeans%potential > 0._RK
        assertion = assertion .and. Kmeans%Err%stat /= 1_IK .and. Kmeans%Err%stat /= 2_IK
        assertion = assertion .and. all(Kmeans%Membership > 0_IK) .and. all(Kmeans%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np
        assertion = assertion .and. all( Point(:,Kmeans%Index) == TestData%Point )
        do ic = 1, nc
            assertion = assertion .and. all( Kmeans%Membership(Kmeans%CumSumSize(ic-1)+1:Kmeans%CumSumSize(ic)) == Kmeans%Membership(Kmeans%CumSumSize(ic)) )
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership, mask = Kmeans%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership, mask = Kmeans%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq, mask = Kmeans%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getKmeans_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Calling the Kmeans routine repeatedly should not cause any errors.
    !> This test is also used for benchmarking the performances of different implementations of the Kmeans algorithm.
    function test_benchmark_1() result(assertion)
        use Statistics_mod, only: getRandInt
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: ncMax = 3_IK
        integer(IK) , parameter     :: nfailMax = 100_IK
        real(RK)    , parameter     :: relTol = 1.e-8_RK
        logical                     :: assertion

        type(Kmeans_type)           :: Kmeans
        integer                     :: i, np, nc

        assertion = .true.

        do i = 1, 100

            nc = ncMax ! getRandInt(2,ncMax)
            np = TestData%np ! getRandInt(100,TestData%np)
            Kmeans = Kmeans_type( nd = TestData%nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = 1_IK & ! LCOV_EXCL_LINE
                                , Point = TestData%Point(1:TestData%nd,1:np)& ! LCOV_EXCL_LINE
                                , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                                , relTol = relTol & ! LCOV_EXCL_LINE
                                )

            assertion = assertion .and. .not. Kmeans%Err%occurred
            assertion = assertion .and. Kmeans%Err%stat /= 1_IK
            assertion = assertion .and. Kmeans%Err%stat /= 2_IK
            assertion = assertion .and. Kmeans%potential > 0._RK
            assertion = assertion .and. all(Kmeans%Membership(1:np) > 0_IK) .and. all(Kmeans%Membership(1:np) < nc + 1)
            assertion = assertion .and. all(Kmeans%MinDistanceSq(1:np) > 0_IK)
            assertion = assertion .and. all(Kmeans%Size > 0_IK) .and. sum(Kmeans%Size)==TestData%np

            if (Test%isVerboseMode .and. .not. assertion) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Size < 1          =", pack(Kmeans%Size, mask = Kmeans%Size < 1_IK)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership < 1    =", pack(Kmeans%Membership(1:np), mask = Kmeans%Membership(1:np) < 1_IK)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Membership > nc   =", pack(Kmeans%Membership(1:np), mask = Kmeans%Membership(1:np) > nc)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%MinDistanceSq < 0 =", pack(Kmeans%MinDistanceSq(1:np), mask = Kmeans%MinDistanceSq(1:np) < 0._RK)
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%occurred      =", Kmeans%Err%occurred
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%potential         =", Kmeans%potential
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Err%stat          =", Kmeans%Err%stat
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%niter             =", Kmeans%niter
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%nfail             =", Kmeans%nfail
                write(Test%outputUnit,"(*(g0.15,:,' '))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_benchmark_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Kmeans_mod ! LCOV_EXCL_LINE