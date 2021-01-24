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
#else
    use Partition_mod
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

    subroutine readTestData(TestData)
        implicit none
        class(TestData_type), intent(inout) :: TestData
        integer(IK) :: ip
!        integer(IK) :: fileUnit
!        open( file = Test%inDir//"/Test_Partition_mod@points.txt" & ! LCOV_EXCL_LINE
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
        integer                     :: fileUnit
        logical                     :: isInside

        assertion = .true.

        Partition = Partition_type  ( Point = TestData%Point & ! LCOV_EXCL_LINE
                                    , nd = TestData%nd & ! LCOV_EXCL_LINE
                                    , np = TestData%np & ! LCOV_EXCL_LINE
                                    !, nemax = TestData%nemax & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedRecursion = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedFailure = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionStabilizationRequested  = .true. & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                    !, mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
#endif
                                    , logTightness = log(1._RK) & ! LCOV_EXCL_LINE
                                    !, inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                    !, parLogVol = sum(log(TestData%DomainSize)) & ! LCOV_EXCL_LINE
                                    , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                    )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_"//MODULE_NAME(2:)//"@test_runPartition_1.Partition."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )

        call Partition%write(fileUnit, TestData%nd, TestData%np, TestData%Point)

        assertion = assertion .and. .not. Partition%Err%occurred
        assertion = assertion .and. Partition%Err%stat /= 1_IK
        assertion = assertion .and. Partition%Err%stat /= 2_IK
        assertion = assertion .and. all(Partition%Size > 0_IK)
        assertion = assertion .and. all(Partition%Membership > 0_IK) .and. all(Partition%Membership <= Partition%neopt)

        do ip = 1, Partition%np
            ic = Partition%Membership(ip)
            NormedPoint = TestData%Point(:,ip) - Partition%Center(:,ic)
            mahalSq = dot_product(NormedPoint,matmul(Partition%InvCovMat(:,:,ic),NormedPoint))
            isInside = mahalSq - 1._RK <= 1.e-6_RK
            if (.not. isInside) then
                if (Test%isVerboseMode) write(Test%outputUnit,"(*(g0.15,:,' '))") "Point not inside!, mahalSq = ", mahalSq
                !assertion = assertion .and. isInside
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
        use Statistics_mod, only: ClusteredPoint_type
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
        integer                     :: fileUnit
        logical                     :: isInside
        integer                     :: ip, ic

        assertion = .true.

        RandomSeed = RandomSeed_type(imageID = Test%Image%id, inputSeed = 1234_IK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate clustered points
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call ClusteredPoint%get ( nd = 2_IK & ! LCOV_EXCL_LINE
                                , etamin = 0.5_RK & ! LCOV_EXCL_LINE
                                , etamax = 10.0_RK & ! LCOV_EXCL_LINE
                                , sizeMin = 100 & ! LCOV_EXCL_LINE
                                , sizeMax = 100 & ! LCOV_EXCL_LINE
                                , centerMin = -2._RK & ! LCOV_EXCL_LINE
                                , centerMax = +2._RK & ! LCOV_EXCL_LINE
                                , nc = 5_IK & ! LCOV_EXCL_LINE
                                !, Size = [50, 1000, 500, 2000, 3] & ! LCOV_EXCL_LINE
                                !, Eta = [1._RK, 2._RK, 0.5_RK, 0.05_RK, 1.5_RK] & ! LCOV_EXCL_LINE
                                !, dist = "uniform" & ! LCOV_EXCL_LINE
                                )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_"//MODULE_NAME(2:)//"@test_runPartition_2.ClusteredPoint."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
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

        RandomSeed = RandomSeed_type(imageID = Test%Image%id, inputSeed = 12243_IK)

        Partition = Partition_type  ( Point = ClusteredPoint%Point & ! LCOV_EXCL_LINE
                                    , nd = ClusteredPoint%nd & ! LCOV_EXCL_LINE
                                    , np = ClusteredPoint%np & ! LCOV_EXCL_LINE
                                    !, nemax = 2_IK & ! LCOV_EXCL_LINE
                                    !, nemax = ClusteredPoint%nemax & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedFailure = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedRecursion = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionMaxAllowedKmeansFailure = 10000 & ! LCOV_EXCL_LINE
                                    !, partitionStabilizationRequested = .true. & ! LCOV_EXCL_LINE
                                    , partitionOptimizationRequested = .true. & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                    , mahalSqWeightExponent = 1._RK & ! LCOV_EXCL_LINE
#endif
                                    , inclusionFraction = 0._RK & ! LCOV_EXCL_LINE
                                    , logTightness = log(1.5_RK) & ! LCOV_EXCL_LINE
                                    !, parLogVol = sum(log(ClusteredPoint%DomainSize)) & ! LCOV_EXCL_LINE
                                    , trimEnabled = .true. & ! LCOV_EXCL_LINE
                                    )

        ! write data to output for further investigation

        open( file = Test%outDir//"/Test_"//MODULE_NAME(2:)//"@test_runPartition_2.Partition."//num2str(Test%Image%id)//".txt" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            , newunit = fileUnit & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        call Partition%write(fileUnit, ClusteredPoint%nd, ClusteredPoint%np, ClusteredPoint%Point)
        close(fileUnit)

!write(*,"(40(g0,:,','))") 
!write(*,"(40(g0,:,','))") 
!write(*,"(40(g0,:,','))") 
!write(*,"(40(g0,:,','))") "size(PartitionMembership), np", size(Partition%Membership), Partition%np
!write(*,"(40(g0,:,','))") Partition%Membership
!write(*,"(40(g0,:,','))") 
!write(*,"(40(g0,:,','))") 
!write(*,"(40(g0,:,','))") 

        assertion = assertion .and. .not. Partition%Err%occurred
        assertion = assertion .and. Partition%Err%stat /= 1_IK
        assertion = assertion .and. Partition%Err%stat /= 2_IK
        assertion = assertion .and. all(Partition%Size > 0_IK)
        assertion = assertion .and. all(Partition%Membership > 0_IK) .and. all(Partition%Membership < Partition%neopt + 1)

        do ip = 1, Partition%np
            ic = Partition%Membership(ip)
            NormedPoint = ClusteredPoint%Point(:,ip) - Partition%Center(:,ic)
            mahalSq = dot_product(NormedPoint,matmul(Partition%InvCovMat(:,:,ic),NormedPoint))
            isInside = mahalSq - 1._RK <= 1.e-6_RK
            if (.not. isInside) then
                if (Test%isVerboseMode) write(Test%outputUnit,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                !assertion = assertion .and. isInside
                exit
            end if
        end do

        !assertion = assertion .and. isInside

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

    end function test_runPartition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%