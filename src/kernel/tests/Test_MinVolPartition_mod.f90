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

!>  \brief This module contains tests of the module [MinVolPartition_mod](@ref MinVolPartition_mod).
!>  \author Amir Shahmoradi

module Test_MinVolPartition_mod

    use Test_mod, only: Test_type
    use MinVolPartition_mod
    implicit none

    private
    public :: test_MinVolPartition

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

    subroutine test_MinVolPartition()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call TestData%read()
        call Test%run(test_runMinVolPartition_1, "test_runMinVolPartition_1")
        call Test%run(test_runMinVolPartition_2, "test_runMinVolPartition_2")
        !call Test%run(test_runMinVolPartition_3, "test_runMinVolPartition_3")
        !call Test%run(test_runMinVolPartition_4, "test_runMinVolPartition_4")
        !call Test%run(test_benchmark_1, "test_benchmark_1")
        !call Test%run(test_optimizeMinVolPartition_1, "test_optimizeMinVolPartition_1")
        call Test%finalize()
    end subroutine test_MinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine readTestData(TestData)
        implicit none
        class(TestData_type), intent(inout) :: TestData
        integer(IK) :: ip
!        integer(IK) :: fileUnit
!        open( file = Test%inDir//"/Test_MinVolPartition_mod@points.txt" & ! LCOV_EXCL_LINE
!            , newunit = fileUnit & ! LCOV_EXCL_LINE
!            , status = "old" & ! LCOV_EXCL_LINE
!#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
!            , SHARED & ! LCOV_EXCL_LINE
!#endif
!            )
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
!        read(fileUnit,*)
!        do ip = 1, TestData%np
!            read(fileUnit,*) TestData%Point(1:TestData%nd,ip)
!        end do
!        close(fileUnit)
    end subroutine readTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runMinVolPartition()` by passing a fixed initial set of cluster centers to the MinVolPartition constructor.
    function test_runMinVolPartition_1() result(assertion)
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion

        type(MinVolPartition_type)  :: MinVolPartition
        integer                     :: fileUnit

        assertion = .true.

        MinVolPartition = MinVolPartition_type  ( Point = TestData%Point & ! LCOV_EXCL_LINE
                                                , nd = TestData%nd & ! LCOV_EXCL_LINE
                                                , np = TestData%np & ! LCOV_EXCL_LINE
                                                !, nemax = TestData%nemax & ! LCOV_EXCL_LINE
                                                !, maxKvolumeLoopRecursion = 10000 & ! LCOV_EXCL_LINE
                                                !, maxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                                !, maxAllowedMinVolFailure = 10000 & ! LCOV_EXCL_LINE
                                                !, stabilizationRequested  = .true. & ! LCOV_EXCL_LINE
                                                !, mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
                                                , tightness  = 1._RK & ! LCOV_EXCL_LINE
                                                !, inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                                !, parLogVol = sum(log(TestData%DomainSize)) & ! LCOV_EXCL_LINE
                                                , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                                )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_MinVolPartition_mod@test_runMinVolPartition_1."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )

        call MinVolPartition%write(fileUnit, TestData%nd, TestData%np, TestData%Point)

        assertion = assertion .and. .not. MinVolPartition%Err%occurred
        assertion = assertion .and. MinVolPartition%Err%stat /= 1_IK
        assertion = assertion .and. MinVolPartition%Err%stat /= 2_IK
        assertion = assertion .and. all(MinVolPartition%Size > 0_IK)
        assertion = assertion .and. all(MinVolPartition%Membership > 0_IK) .and. all(MinVolPartition%Membership < MinVolPartition%neopt + 1)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%neopt                =", MinVolPartition%neopt
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Size < 1             =", pack(MinVolPartition%Size, mask = MinVolPartition%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Membership < 1       =", pack(MinVolPartition%Membership, mask = MinVolPartition%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Membership > neopt   =", pack(MinVolPartition%Membership, mask = MinVolPartition%Membership > MinVolPartition%neopt)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%occurred         =", MinVolPartition%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%stat             =", MinVolPartition%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%msg              =", MinVolPartition%Err%msg
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runMinVolPartition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> test `runMinVolPartition()` by passing a fixed initial set of cluster centers to the MinVolPartition constructor.
    function test_runMinVolPartition_2() result(assertion)
        use Statistics_mod, only: ClusteredPoint_type
        use RandomSeed_mod, only: RandomSeed_type
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        type(MinVolPartition_type)  :: MinVolPartition
        type(ClusteredPoint_type)   :: ClusteredPoint
        type(RandomSeed_type)       :: RandomSeed
        integer                     :: fileUnit

        assertion = .true.

        !RandomSeed = RandomSeed_type(imageID = Test%Image%id, inputSeed = 1234_IK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call ClusteredPoint%get ( & ! LCOV_EXCL_LINE
                                 nd = 2_IK & ! LCOV_EXCL_LINE
                                !, nc = 5_IK, & ! LCOV_EXCL_LINE
                                !, Size = [50, 1000, 500, 2000, 3] & ! LCOV_EXCL_LINE
                                !, Eta = [1._RK, 2._RK, 0.5_RK, 0.05_RK, 1.5_RK] & ! LCOV_EXCL_LINE
                                , dist = "uniform" & ! LCOV_EXCL_LINE
                                )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_MinVolPartition_mod@test_runMinVolPartition_2.ClusteredPoint."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        call ClusteredPoint%write(fileUnit)
        close(fileUnit)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Partition the clustered points clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        MinVolPartition = MinVolPartition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                                , nd = ClusteredPoint%nd & ! LCOV_EXCL_LINE
                                                , np = ClusteredPoint%np & ! LCOV_EXCL_LINE
                                                !, nemax = ClusteredPoint%nemax & ! LCOV_EXCL_LINE
                                                !, maxKvolumeLoopRecursion = 10000 & ! LCOV_EXCL_LINE
                                                !, maxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                                !, maxAllowedMinVolFailure = 10000 & ! LCOV_EXCL_LINE
                                                !, stabilizationRequested  = .true. & ! LCOV_EXCL_LINE
                                                !, mahalSqWeightExponent = 0._RK & ! LCOV_EXCL_LINE
                                                !, tightness  = 1._RK & ! LCOV_EXCL_LINE
                                                !, inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                                !, parLogVol = sum(log(ClusteredPoint%DomainSize)) & ! LCOV_EXCL_LINE
                                                , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                                )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_MinVolPartition_mod@test_runMinVolPartition_2.MinVolPartition."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        call MinVolPartition%write(fileUnit, ClusteredPoint%nd, ClusteredPoint%np, ClusteredPoint%Point)
        close(fileUnit)

        assertion = assertion .and. .not. MinVolPartition%Err%occurred
        assertion = assertion .and. MinVolPartition%Err%stat /= 1_IK
        assertion = assertion .and. MinVolPartition%Err%stat /= 2_IK
        assertion = assertion .and. all(MinVolPartition%Size > 0_IK)
        assertion = assertion .and. all(MinVolPartition%Membership > 0_IK) .and. all(MinVolPartition%Membership < MinVolPartition%neopt + 1)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%neopt                =", MinVolPartition%neopt
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Size < 1             =", pack(MinVolPartition%Size, mask = MinVolPartition%Size < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Membership < 1       =", pack(MinVolPartition%Membership, mask = MinVolPartition%Membership < 1_IK)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Membership > neopt   =", pack(MinVolPartition%Membership, mask = MinVolPartition%Membership > MinVolPartition%neopt)
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%occurred         =", MinVolPartition%Err%occurred
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%stat             =", MinVolPartition%Err%stat
            write(Test%outputUnit,"(*(g0.15,:,' '))") "MinVolPartition%Err%msg              =", MinVolPartition%Err%msg
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_runMinVolPartition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_MinVolPartition_mod ! LCOV_EXCL_LINE