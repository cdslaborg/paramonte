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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Kmeans()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_runKmeans_1, "test_runKmeans_1")
        call Test%run(test_runKmeans_2, "test_runKmeans_2")
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
        type(TestData_type)         :: TestData
        type(Kmeans_type)           :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call TestData%read()

        Kmeans = Kmeans_type( nc = nc & ! LCOV_EXCL_LINE
                            , nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            )
        if (.not. assertion) return

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_Kmeans_mod@test_runKmeans_1@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_Kmeans_mod@test_runKmeans_1@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Center(1:TestData%nd,i), Kmeans%Size(i)
        end do
        close(fileUnit)

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
        type(TestData_type)         :: TestData
        type(Kmeans_type)           :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call TestData%read()

        InitCenter = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])

        Kmeans = Kmeans_type( nc = nc & ! LCOV_EXCL_LINE
                            , nd = TestData%nd & ! LCOV_EXCL_LINE
                            , np = TestData%np & ! LCOV_EXCL_LINE
                            , Point = TestData%Point & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )
        if (.not. assertion) return

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_Kmeans_mod@test_runKmeans_2@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_Kmeans_mod@test_runKmeans_2@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Center(1:TestData%nd,i), Kmeans%Size(i)
        end do
        close(fileUnit)

    end function test_runKmeans_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Kmeans_mod