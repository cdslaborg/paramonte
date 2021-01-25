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

!>  \brief
!>  This module contains tests of the module [KmeansOOP_mod](@ref kmeansoop_mod).
!>  The OOP version of Kmeans clustering algorithm happens to be 100-300% slower than
!>  the functional version as implemented in the module [Kmeans_mod](@ref kmeans_mod).
!>  \author Amir Shahmoradi

module Test_KmeansOOP_mod

    use Test_mod, only: Test_type
    use KmeansOOP_mod
    implicit none

    private
    public :: test_KmeansOOP

    type(Test_type) :: Test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 1000
        real(RK), allocatable   :: Point(:,:)
    contains
        procedure, pass :: read => readTestData
    end type TestData_type

    type(TestData_type)         :: TestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_KmeansOOP()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call TestData%read()
        call Test%run(test_runKmeans_1, "test_runKmeans_1")
        call Test%run(test_runKmeans_2, "test_runKmeans_2")
        call Test%run(test_runKmeans_3, "test_runKmeans_3")
        call Test%run(test_runKmeans_4, "test_runKmeans_4")
        call Test%run(test_runKmeans_5, "test_runKmeans_5")
        call Test%run(test_benchmark_1, "test_benchmark_1")
        call Test%finalize()
    end subroutine test_KmeansOOP

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
        integer(IK), parameter      :: nt = 1_IK
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call Kmeans%run ( nd = TestData%nd          & ! LCOV_EXCL_LINE
                        , np = TestData%np          & ! LCOV_EXCL_LINE
                        , nc = nc                   & ! LCOV_EXCL_LINE
                        , nt = nt                   & ! LCOV_EXCL_LINE
                        , Point = TestData%Point    & ! LCOV_EXCL_LINE
                        )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_1@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_1@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Try(Kmeans%itopt)%Center(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Size(i)
        end do
        close(fileUnit)

        assertion = assertion .and. .not. Kmeans%Try(Kmeans%itopt)%Err%occurred
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Membership > 0_IK) .and. all(Kmeans%Try(Kmeans%itopt)%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Size > 0_IK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Size < 1          =", pack(Kmeans%Try(Kmeans%itopt)%Size, mask = Kmeans%Try(Kmeans%itopt)%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership < 1    =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership > nc   =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0 =", pack(Kmeans%Try(Kmeans%itopt)%MinDistanceSq, mask = Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%occurred      =", Kmeans%Try(Kmeans%itopt)%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%potential         =", Kmeans%Try(Kmeans%itopt)%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%stat          =", Kmeans%Try(Kmeans%itopt)%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%niter             =", Kmeans%Try(Kmeans%itopt)%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%nzsci             =", Kmeans%Try(Kmeans%itopt)%nzsci
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
        integer(IK), parameter      :: nt = 1_IK
        real(RK)    , allocatable   :: InitCenter(:,:,:)
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        allocate(InitCenter(TestData%nd,nc,nt))
        InitCenter(:,:,1) = reshape([4.7_RK, 4.7_RK, 6.4_RK, 6.1_RK, 9.5_RK, 8.6_RK], shape = [TestData%nd,nc])
        do i = 2, nt
            InitCenter(:,:,i) = InitCenter(:,:,1)
        end do

        call Kmeans%run ( nd = TestData%nd          & ! LCOV_EXCL_LINE
                        , np = TestData%np          & ! LCOV_EXCL_LINE
                        , nc = nc                   & ! LCOV_EXCL_LINE
                        , nt = nt                   & ! LCOV_EXCL_LINE
                        , Point = TestData%Point    & ! LCOV_EXCL_LINE
                        , InitCenter = InitCenter   & ! LCOV_EXCL_LINE
                        )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_2@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_2@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Try(Kmeans%itopt)%Center(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Size(i)
        end do
        close(fileUnit)

        assertion = assertion .and. .not. Kmeans%Try(Kmeans%itopt)%Err%occurred
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Membership > 0_IK) .and. all(Kmeans%Try(Kmeans%itopt)%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Size > 0_IK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Size < 1          =", pack(Kmeans%Try(Kmeans%itopt)%Size, mask = Kmeans%Try(Kmeans%itopt)%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership < 1    =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership > nc   =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0 =", pack(Kmeans%Try(Kmeans%itopt)%MinDistanceSq, mask = Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%occurred      =", Kmeans%Try(Kmeans%itopt)%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%potential         =", Kmeans%Try(Kmeans%itopt)%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%stat          =", Kmeans%Try(Kmeans%itopt)%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%niter             =", Kmeans%Try(Kmeans%itopt)%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%nzsci             =", Kmeans%Try(Kmeans%itopt)%nzsci
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
        integer(IK), parameter      :: nt = 1_IK
        integer(IK) , parameter     :: niterMax = 1_IK
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call Kmeans%run ( nd = TestData%nd          & ! LCOV_EXCL_LINE
                        , np = TestData%np          & ! LCOV_EXCL_LINE
                        , nc = nc                   & ! LCOV_EXCL_LINE
                        , nt = nt                   & ! LCOV_EXCL_LINE
                        , niterMax = niterMax       & ! LCOV_EXCL_LINE
                        , Point = TestData%Point    & ! LCOV_EXCL_LINE
                        )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_3@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Try(1)%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_3@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Try(1)%Center(1:TestData%nd,i), Kmeans%Try(1)%Size(i)
        end do
        close(fileUnit)

        assertion = assertion .and. Kmeans%failed
        assertion = assertion .and. Kmeans%itopt == 0_IK
        assertion = assertion .and. Kmeans%Try(1)%Err%occurred
        assertion = assertion .and. Kmeans%Try(1)%Err%stat == 1_IK
        assertion = assertion .and. Kmeans%Try(1)%niter == niterMax
        assertion = assertion .and. Kmeans%Try(1)%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Try(1)%Membership > 0_IK) .and. all(Kmeans%Try(1)%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%Try(1)%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Try(1)%Size > 0_IK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Size < 1           =", pack(Kmeans%Try(1)%Size, mask = Kmeans%Try(1)%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Membership < 1     =", pack(Kmeans%Try(1)%Membership, mask = Kmeans%Try(1)%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Membership > nc    =", pack(Kmeans%Try(1)%Membership, mask = Kmeans%Try(1)%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%MinDistanceSq < 0  =", pack(Kmeans%Try(1)%MinDistanceSq, mask = Kmeans%Try(1)%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Err%occurred       =", Kmeans%Try(1)%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%potential          =", Kmeans%Try(1)%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Err%stat           =", Kmeans%Try(1)%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%niter              =", Kmeans%Try(1)%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%nzsci              =", Kmeans%Try(1)%nzsci
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%failed                    =", Kmeans%failed
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%itopt                     =", Kmeans%itopt
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The function `runKmeans()` must function properly for reasonable optional input values of `nzsciMax` and `reltolSq`.
    function test_runKmeans_4() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK), parameter      :: nt = 1_IK
        integer(IK) , parameter     :: nzsciMax = 100_IK
        real(RK)    , parameter     :: reltol = 1.e-8_RK
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call Kmeans%run ( nd = TestData%nd          & ! LCOV_EXCL_LINE
                        , np = TestData%np          & ! LCOV_EXCL_LINE
                        , nc = nc                   & ! LCOV_EXCL_LINE
                        , nt = nt                   & ! LCOV_EXCL_LINE
                        , reltol = reltol           & ! LCOV_EXCL_LINE
                        , nzsciMax = nzsciMax       & ! LCOV_EXCL_LINE
                        , Point = TestData%Point    & ! LCOV_EXCL_LINE
                        )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_4@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_4@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Try(1)%Center(1:TestData%nd,i), Kmeans%Try(1)%Size(i)
        end do
        close(fileUnit)

        assertion = assertion .and. .not. Kmeans%failed
        assertion = assertion .and. .not. Kmeans%Try(1)%Err%occurred
        assertion = assertion .and. Kmeans%Try(1)%Err%stat /= 1_IK
        assertion = assertion .and. Kmeans%Try(1)%Err%stat /= 2_IK
        assertion = assertion .and. Kmeans%Try(1)%potential > 0._RK
        assertion = assertion .and. all(Kmeans%Try(1)%Membership > 0_IK) .and. all(Kmeans%Try(1)%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%Try(1)%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Try(1)%Size > 0_IK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Size < 1           =", pack(Kmeans%Try(1)%Size, mask = Kmeans%Try(1)%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Membership < 1     =", pack(Kmeans%Try(1)%Membership, mask = Kmeans%Try(1)%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Membership > nc    =", pack(Kmeans%Try(1)%Membership, mask = Kmeans%Try(1)%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%MinDistanceSq < 0  =", pack(Kmeans%Try(1)%MinDistanceSq, mask = Kmeans%Try(1)%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Err%occurred       =", Kmeans%Try(1)%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%potential          =", Kmeans%Try(1)%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%Err%stat           =", Kmeans%Try(1)%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%niter              =", Kmeans%Try(1)%niter
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(1)%nzsci              =", Kmeans%Try(1)%nzsci
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%failed                    =", Kmeans%failed
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%itopt                     =", Kmeans%itopt
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runKmeans()` by passing a number of tries to find the more optimal Kmeans clustering.
    function test_runKmeans_5() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: nc = 3_IK
        integer(IK) , parameter     :: nt = 2_IK + nint(log(real(nc)))
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, fileUnit

        assertion = .true.

        call Kmeans%run ( nd = TestData%nd          & ! LCOV_EXCL_LINE
                        , np = TestData%np          & ! LCOV_EXCL_LINE
                        , nc = nc                   & ! LCOV_EXCL_LINE
                        , nt = nt                   & ! LCOV_EXCL_LINE
                        , Point = TestData%Point    & ! LCOV_EXCL_LINE
                        )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_5@points."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "membership"
        do i = 1, TestData%np
            write(fileUnit,"(*(g0.15,:,','))") TestData%Point(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Membership(i)
        end do
        close(fileUnit)

        open( file = Test%outDir//"/Test_KmeansOOP_mod@test_runKmeans_5@centers."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        write(fileUnit,"(*(g0.15,:,' '))") "x", "y", "size"
        do i = 1, nc
            write(fileUnit,"(*(g0.15,:,','))") Kmeans%Try(Kmeans%itopt)%Center(1:TestData%nd,i), Kmeans%Try(Kmeans%itopt)%Size(i)
        end do
        close(fileUnit)

        assertion = assertion .and. .not. Kmeans%Try(Kmeans%itopt)%Err%occurred
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%potential > 0._RK
        assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 1_IK .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 2_IK
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Membership > 0_IK) .and. all(Kmeans%Try(Kmeans%itopt)%Membership < nc + 1)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%MinDistanceSq > 0_IK)
        assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Size > 0_IK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Size < 1          =", pack(Kmeans%Try(Kmeans%itopt)%Size, mask = Kmeans%Try(Kmeans%itopt)%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership < 1    =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership > nc   =", pack(Kmeans%Try(Kmeans%itopt)%Membership, mask = Kmeans%Try(Kmeans%itopt)%Membership > nc)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0 =", pack(Kmeans%Try(Kmeans%itopt)%MinDistanceSq, mask = Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0._RK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%occurred      =", Kmeans%Try(Kmeans%itopt)%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%potential         =", Kmeans%Try(Kmeans%itopt)%potential
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%stat          =", Kmeans%Try(Kmeans%itopt)%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runKmeans_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Calling the Kmeans routine repeatedly should not cause any errors.
    !> This test is also used for benchmarking the performances of different implementations of the Kmeans algorithm.
    function test_benchmark_1() result(assertion)
        use Statistics_mod, only: getRandInt
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        integer(IK) , parameter     :: ncMax = 3_IK
        integer(IK), parameter      :: ntMax = 1_IK
        logical                     :: assertion
        type(KmeansOOP_type)        :: Kmeans
        integer                     :: i, np, nc

        assertion = .true.

        Kmeans = KmeansOOP_type(TestData%nd, TestData%np, ncMax, ntMax)
        do i = 1, 1000

            nc = ncMax ! getRandInt(2,ncMax)
            np = TestData%np ! getRandInt(100,TestData%np)
            call Kmeans%run ( nd = TestData%nd                          & ! LCOV_EXCL_LINE
                            , np = np                                   & ! LCOV_EXCL_LINE
                            , nc = nc                                   & ! LCOV_EXCL_LINE
                            , nt = ntMax                                & ! LCOV_EXCL_LINE
                            , Point = TestData%Point(1:TestData%nd,1:np)& ! LCOV_EXCL_LINE
                            )

            if (.not. Kmeans%failed) then
                assertion = assertion .and. .not. Kmeans%Try(Kmeans%itopt)%Err%occurred
                assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 1_IK
                assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%Err%stat /= 2_IK
                assertion = assertion .and. Kmeans%Try(Kmeans%itopt)%potential > 0._RK
                assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Membership(1:np) > 0_IK) .and. all(Kmeans%Try(Kmeans%itopt)%Membership(1:np) < nc + 1)
                assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%MinDistanceSq(1:np) > 0_IK)
                assertion = assertion .and. all(Kmeans%Try(Kmeans%itopt)%Size > 0_IK)

                if (Test%isVerboseMode .and. .not. assertion) then
                ! LCOV_EXCL_START
                    write(Test%outputUnit,"(*(g0.15,:,' '))")
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Size < 1            =", pack(Kmeans%Try(Kmeans%itopt)%Size, mask = Kmeans%Try(Kmeans%itopt)%Size < 1_IK)
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership < 1      =", pack(Kmeans%Try(Kmeans%itopt)%Membership(1:np), mask = Kmeans%Try(Kmeans%itopt)%Membership(1:np) < 1_IK)
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Membership > nc     =", pack(Kmeans%Try(Kmeans%itopt)%Membership(1:np), mask = Kmeans%Try(Kmeans%itopt)%Membership(1:np) > nc)
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%MinDistanceSq < 0   =", pack(Kmeans%Try(Kmeans%itopt)%MinDistanceSq(1:np), mask = Kmeans%Try(Kmeans%itopt)%MinDistanceSq(1:np) < 0._RK)
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%occurred        =", Kmeans%Try(Kmeans%itopt)%Err%occurred
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%potential           =", Kmeans%Try(Kmeans%itopt)%potential
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%Err%stat            =", Kmeans%Try(Kmeans%itopt)%Err%stat
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%niter               =", Kmeans%Try(Kmeans%itopt)%niter
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%Try(Kmeans%itopt)%nzsci               =", Kmeans%Try(Kmeans%itopt)%nzsci
                    write(Test%outputUnit,"(*(g0.15,:,' '))") "Kmeans%itopt                                 =", Kmeans%itopt
                    write(Test%outputUnit,"(*(g0.15,:,' '))")
                    !error stop
                end if
                ! LCOV_EXCL_STOP
            end if

        end do

    end function test_benchmark_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_KmeansOOP_mod ! LCOV_EXCL_LINE