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

!>  \brief This module contains tests of the module [KnnMaxDen_mod](@ref Knnmaxden_mod).
!>  \author Amir Shahmoradi

module Test_Knn_mod

!>  \brief This file implements the bodies of the Knn module.
!>  \author Amir Shahmoradi

    use Test_mod, only: Test_type

    use JaggedArray_mod, only: IV => IntVec_type, RV => RealVec_type
    use Knn_mod

    implicit none

    private
    public :: test_Knn

    type(Test_type) :: Test

    type :: TestData_type
        integer(IK)             :: nd = 2
        integer(IK)             :: np = 10
        integer(IK)             :: HubMinDistEdgeIndex(10) = [5, 8, 7, 7, 10, 9, 3, 2, 3, 5]
        real(RK)                :: HubMinDistEdgeLenSq(10) =    [ 0.693595837531895E-1_RK &
                                                                , 0.411631653715387E-2_RK &
                                                                , 0.324022691025059E-2_RK &
                                                                , 0.122396412769469E-1_RK &
                                                                , 0.106512144573307E-1_RK &
                                                                , 0.442485474967571_RK &
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

    subroutine test_Knn()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        TestData = TestData_type()
        call Test%run(test_hubify_1, "test_hubify_1")
        call Test%run(test_getDistSq_1, "test_getDistSq_1")
        call Test%run(test_getLogDensity_1, "test_getLogDensity_1")
        call Test%finalize()
    end subroutine test_Knn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function contructTestData() result(TestData)
        implicit none
        type(TestData_type) :: TestData
        allocate(TestData%HubEdgeIndex(7))
        TestData%HubEdgeIndex(1)%Vector = [7, 9]
        TestData%HubEdgeIndex(2)%Vector = [3, 4]
        TestData%HubEdgeIndex(3)%Vector = [1,10]
        TestData%HubEdgeIndex(4)%Vector = [8]
        TestData%HubEdgeIndex(5)%Vector = [6]
        TestData%HubEdgeIndex(6)%Vector = [5]
        TestData%HubEdgeIndex(7)%Vector = [2]
        allocate(TestData%HubEdgeLenSq(7))
        TestData%HubEdgeLenSq(1)%Vector = [0.324022691025059E-2_RK, 0.375960647686508E-1_RK]
        TestData%HubEdgeLenSq(2)%Vector = [0.324022691025059E-2_RK, 0.122396412769469E-1_RK]
        TestData%HubEdgeLenSq(3)%Vector = [0.693595837531895E-1_RK, 0.106512144573307E-1_RK]
        TestData%HubEdgeLenSq(4)%Vector = [0.411631653715387E-2_RK]
        TestData%HubEdgeLenSq(5)%Vector = [0.442485474967571_RK]
        TestData%HubEdgeLenSq(6)%Vector = [0.106512144573307E-1_RK]
        TestData%HubEdgeLenSq(7)%Vector = [0.411631653715387E-2_RK]
    end function contructTestData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getDistSq()`.
    function test_getDistSq_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK), allocatable       :: Dist(:,:)
        real(RK), allocatable       :: Diff(:,:)
        real(RK), parameter         :: tolerance = 1.e-10_RK
        integer                     :: ip, jp

        assertion = .true.

        Dist = sqrt(getDistSq(TestData%nd, TestData%np, TestData%Point))
        Diff = abs( (Dist - TestData%Dist) )

        assertion = all( Diff < tolerance )

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Dist, Dist_ref, Diff"
            do jp = 1, TestData%np
                do ip = 1, TestData%np
                    write(Test%outputUnit,"(*(g0.15,:,' '))") Dist(ip,jp), TestData%Dist(ip,jp), Diff(ip,jp)
                end do
            end do
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDistSq_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `getMeanMinDist()`.
    function test_getLogDensity_1() result(assertion)

        use ClusteredPoint_mod, only: ClusteredPoint_type
        use Constants_mod, only: IK, RK
        use Math_mod, only: getLogVolUnitBall

        implicit none
        logical                     :: assertion

        type(ClusteredPoint_type)   :: ClusteredPoint
        real(RK)                    :: logDensityFromMeanMinDist
        real(RK)                    :: logVolUnitBall
        real(RK)                    :: logMeanMinDist
        integer                     :: i, ih, nh

        character(:), allocatable   :: dist
        integer(IK)                 :: rngseed, nd, nc, sizeMin, sizeMax
        real(RK)                    :: etamin, etamax, centerMin, centerMax
        namelist /specData/ rngseed, nd, nc, sizeMin, sizeMax, etamin, etamax, centerMin, centerMax, dist

        assertion = .true.

        rngseed = -huge(rngseed)
        allocate(character(63) :: dist)
        open(newunit = Test%File%unit, file = Test%inDir//"/Test_Knn_mod@test_getLogDensity_1.nml", status = "old")
        read(Test%File%unit, nml = specData)
        dist = trim(adjustl(dist))
        close(Test%File%unit)

        logVolUnitBall = getLogVolUnitBall(nd)

!write(*,*) "logVolUnitBall", logVolUnitBall

        Test%File = Test%openFile() ! label = "LogVolEstimate")

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
                                    )

            ! write data to output for further investigation

            !Test%File = Test%openFile(label = "ClusteredPoint")
            !call ClusteredPoint%write(Test%File%unit)
            !close(Test%File%unit)

            !logMeanMinDist = 0._RK
            !nh = ClusteredPoint%Hub%nh / 2
            !do ih = 1, nh
            !    logMeanMinDist = logMeanMinDist + sum(sqrt(ClusteredPoint%Hub%EdgeLenSq(ih)%Vector(1:1)))
            !    if (size(ClusteredPoint%Hub%EdgeLenSq(ih)%Vector) /= ClusteredPoint%Hub%EdgeCount(ih)) error stop
            !end do
            !logMeanMinDist = log( logMeanMinDist / nh ) ! / ClusteredPoint%Hub%EdgeCount(1:2)

            !logMeanMinDist = log( sum(sqrt(ClusteredPoint%Hub%EdgeLenSq(1)%Vector)) / ClusteredPoint%Hub%EdgeCount(1) )
            !if (ClusteredPoint%Hub%EdgeCount(1) /= size((ClusteredPoint%Hub%EdgeLenSq(1)%Vector))) error stop

            logMeanMinDist = log( sqrt(minval(ClusteredPoint%Hub%EdgeLenSq(1)%Vector)) )

            logDensityFromMeanMinDist = getLogDensity(nd = nd, logMeanMinDist = logMeanMinDist, logVolUnitBall = logVolUnitBall)

            write(Test%File%unit, "(*(g0,:,','))") ClusteredPoint%LogDensity, logDensityFromMeanMinDist

        end do

        close(Test%File%unit)

    end function test_getLogDensity_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Test `hubify()`.
    function test_hubify_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                     :: assertion
        real(RK), parameter         :: tolerance = 1.e-10_RK
        type(Hub_type)              :: Hub
        integer                     :: i

        assertion = .true.

        Hub = Hub_type  ( np = TestData%np & ! LCOV_EXCL_LINE
                        , DistSq = getDistSq(nd = TestData%nd, np = TestData%np, Point = TestData%Point) & ! LCOV_EXCL_LINE
                        )

        assertion = .not. Hub%Err%occurred
        if (.not. assertion) return

        assertion = assertion .and. all( Hub%MinDistEdge%Index == TestData%HubMinDistEdgeIndex )
        if (Test%isVerboseMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%MinDistEdge%Index        ", Hub%MinDistEdge%Index
            write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubMinDistEdgeIndex ", TestData%HubMinDistEdgeIndex
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( abs((Hub%MinDistEdge%LenSq - TestData%HubMinDistEdgeLenSq)/TestData%HubMinDistEdgeLenSq) < tolerance )
        if (Test%isVerboseMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%MinDistEdge%LenSq        ", Hub%MinDistEdge%LenSq
            write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubMinDistEdgeLenSq ", TestData%HubMinDistEdgeLenSq
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( Hub%NodeIndex == TestData%HubNodeIndex )
        if (Test%isVerboseMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%NodeIndex        ", Hub%NodeIndex
            write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubNodeIndex", TestData%HubNodeIndex
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. all( Hub%EdgeCount == TestData%HubEdgeCount )
        if (Test%isVerboseMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%EdgeCount        ", Hub%EdgeCount
            write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubEdgeCount", TestData%HubEdgeCount
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            ! LCOV_EXCL_STOP
        end if
        if (.not. assertion) return ! LCOV_EXCL_LINE

        do i = 1, Hub%nh
            assertion = assertion .and. all(Hub%EdgeIndex(i)%Vector == TestData%HubEdgeIndex(i)%Vector)
            if (Test%isVerboseMode .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%EdgeIndex(i)%Vector          ", Hub%EdgeIndex(i)%Vector
                write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubEdgeIndex(i)%Vector  ", TestData%HubEdgeIndex(i)%Vector
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
            if (.not. assertion) return ! LCOV_EXCL_LINE
        end do

        do i = 1, Hub%nh
            assertion = assertion .and. all( abs((Hub%EdgeLenSq(i)%Vector - TestData%HubEdgeLenSq(i)%Vector)/TestData%HubEdgeLenSq(i)%Vector) < tolerance )
            if (Test%isVerboseMode .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                write(Test%outputUnit,"(*(g0.15,:,' '))") "Hub%EdgeLenSq(i)%Vector          ", Hub%EdgeLenSq(i)%Vector
                write(Test%outputUnit,"(*(g0.15,:,' '))") "TestData%HubEdgeLenSq(i)%Vector  ", TestData%HubEdgeLenSq(i)%Vector
                write(Test%outputUnit,"(*(g0.15,:,' '))")
                ! LCOV_EXCL_STOP
            end if
            if (.not. assertion) return ! LCOV_EXCL_LINE
        end do

    end function test_hubify_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Knn_mod ! LCOV_EXCL_LINE