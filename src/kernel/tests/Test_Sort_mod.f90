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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [Sort_mod](@ref sort_mod).
!>  \author Amir Shahmoradi

module Test_Sort_mod

    use Sort_mod
    use Test_mod, only: Test_type
    use Constants_mod, only: IK, RK
    implicit none

    private
    public :: test_Sort

    type(Test_type) :: Test
    integer, parameter  :: ndata = 50
    integer(IK), parameter :: DataUnsorted_IK(ndata)= &
                                                    [ 1201_IK &
                                                    , 1187_IK &
                                                    , 1188_IK &
                                                    , 1193_IK &
                                                    , 1177_IK &
                                                    , 1153_IK &
                                                    , 1134_IK &
                                                    , 1146_IK &
                                                    , 1172_IK &
                                                    , 1181_IK &
                                                    , 1197_IK &
                                                    , 1172_IK &
                                                    , 1141_IK &
                                                    , 1216_IK &
                                                    , 1158_IK &
                                                    , 1174_IK &
                                                    , 1189_IK &
                                                    , 1211_IK &
                                                    , 1157_IK &
                                                    , 1184_IK &
                                                    , 1177_IK &
                                                    , 1157_IK &
                                                    , 1191_IK &
                                                    , 1176_IK &
                                                    , 1196_IK &
                                                    , 1150_IK &
                                                    , 1185_IK &
                                                    , 1190_IK &
                                                    , 1172_IK &
                                                    , 1161_IK &
                                                    , 1179_IK &
                                                    , 1189_IK &
                                                    , 1136_IK &
                                                    , 1148_IK &
                                                    , 1176_IK &
                                                    , 1142_IK &
                                                    , 1146_IK &
                                                    , 1202_IK &
                                                    , 1156_IK &
                                                    , 1170_IK &
                                                    , 1146_IK &
                                                    , 1170_IK &
                                                    , 1169_IK &
                                                    , 1211_IK &
                                                    , 1168_IK &
                                                    , 1189_IK &
                                                    , 1170_IK &
                                                    , 1162_IK &
                                                    , 1167_IK &
                                                    , 1180_IK ]
    real(RK), parameter :: DataUnsorted_RK(ndata)   = &
                                                    [ 5.28935260000000_RK &
                                                    , 5.50145870000000_RK &
                                                    , 5.89022390000000_RK &
                                                    , 5.06549460000000_RK &
                                                    , 5.62128260000000_RK &
                                                    , 4.49246930000000_RK &
                                                    , 3.54559920000000_RK &
                                                    , 4.17171310000000_RK &
                                                    , 5.34432780000000_RK &
                                                    , 4.30855910000000_RK &
                                                    , 6.12466330000000_RK &
                                                    , 4.45103540000000_RK &
                                                    , 4.08259680000000_RK &
                                                    , 7.64761290000000_RK &
                                                    , 6.53095480000000_RK &
                                                    , 6.07550490000000_RK &
                                                    , 7.32100850000000_RK &
                                                    , 5.82501650000000_RK &
                                                    , 4.19347540000000_RK &
                                                    , 4.89687790000000_RK &
                                                    , 5.61290890000000_RK &
                                                    , 5.70994940000000_RK &
                                                    , 5.00047920000000_RK &
                                                    , 5.47741520000000_RK &
                                                    , 4.99151560000000_RK &
                                                    , 5.08172850000000_RK &
                                                    , 5.98773500000000_RK &
                                                    , 6.97849360000000_RK &
                                                    , 6.91612860000000_RK &
                                                    , 4.90595890000000_RK &
                                                    , 5.71852950000000_RK &
                                                    , 4.12146660000000_RK &
                                                    , 5.51241440000000_RK &
                                                    , 5.26293780000000_RK &
                                                    , 5.14932990000000_RK &
                                                    , 4.14738170000000_RK &
                                                    , 5.55786790000000_RK &
                                                    , 7.08800450000000_RK &
                                                    , 6.08987380000000_RK &
                                                    , 4.73697940000000_RK &
                                                    , 3.80934450000000_RK &
                                                    , 6.03942270000000_RK &
                                                    , 5.96600840000000_RK &
                                                    , 6.06674510000000_RK &
                                                    , 5.84361600000000_RK &
                                                    , 6.19013970000000_RK &
                                                    , 4.43891700000000_RK &
                                                    , 4.45833300000000_RK &
                                                    , 5.47659170000000_RK &
                                                    , 4.65761920000000_RK ]

    integer(IK) , parameter :: DataUnsorted2_IK(ndata) = DataUnsorted_IK
    real(RK)    , parameter :: DataUnsorted2_RK(ndata) = DataUnsorted_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Sort()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_sortArray_1, "test_sortArray_1")
        call Test%run(test_getMedian_RK_1, "test_getMedian_RK_1")
        call Test%run(test_getMedian_RK_2, "test_getMedian_RK_2")
        call Test%run(test_sortAscending_RK_1, "test_sortAscending_RK_1")
        call Test%run(test_sortAscendingWithRooter_IK_1, "test_sortAscendingWithRooter_IK_1")
        call Test%run(test_sortAscendingWithRooter_RK_1, "test_sortAscendingWithRooter_RK_1")
        call Test%finalize()
    end subroutine test_Sort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortArray_1() result(assertion)

        implicit none

        integer(IK)             :: i
        logical                 :: assertion
        real(RK), allocatable   :: DataUnsorted(:)

        assertion = .true.

        DataUnsorted = DataUnsorted_RK
        call sortArray(DataUnsorted)

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted(i-1) <= DataUnsorted(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "DataUnsorted, DataSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted_RK(i), DataUnsorted(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_sortArray_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortAscending_RK_1() result(assertion)

        use Err_mod, only: Err_type
        implicit none

        integer(IK)                 :: i
        logical                     :: assertion
        real(RK)    , allocatable   :: DataUnsorted(:)
        type(Err_type)              :: Err

        DataUnsorted = DataUnsorted_RK
        call sortAscending(np = ndata, Point = DataUnsorted, Err = Err)
        assertion = .not. Err%occurred
        if (.not. assertion) return

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted(i-1) <= DataUnsorted(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "DataUnsorted, DataSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted_RK(i), DataUnsorted(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_sortAscending_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortAscendingWithRooter_IK_1() result(assertion)

        use Err_mod, only: Err_type
        implicit none

        integer(IK)                 :: i
        logical                     :: assertion
        integer(IK) , allocatable   :: DataUnsorted(:)
        integer(IK) , allocatable   :: DataUnsorted2(:)
        type(Err_type)              :: Err

        DataUnsorted    = DataUnsorted_IK
        DataUnsorted2   = DataUnsorted2_IK
        call sortAscendingWithRooter_IK(lenLeader = ndata, Leader = DataUnsorted, Rooter = DataUnsorted2, Err = Err)
        assertion = .not. Err%occurred
        if (.not. assertion) return

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted(i-1) <= DataUnsorted(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "LeaderUnsorted, LeaderSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted_IK(i), DataUnsorted(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted2(i-1) <= DataUnsorted2(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "RooterUnsorted, RooterSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted2_IK(i), DataUnsorted2(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_sortAscendingWithRooter_IK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sortAscendingWithRooter_RK_1() result(assertion)

        use Err_mod, only: Err_type
        implicit none

        integer(IK)                 :: i
        logical                     :: assertion
        real(RK)    , allocatable   :: DataUnsorted(:)
        real(RK)    , allocatable   :: DataUnsorted2(:)
        type(Err_type)              :: Err

        DataUnsorted    = DataUnsorted_RK
        DataUnsorted2   = DataUnsorted2_RK
        call sortAscendingWithRooter_RK(lenLeader = ndata, Leader = DataUnsorted, Rooter = DataUnsorted2, Err = Err)
        assertion = .not. Err%occurred
        if (.not. assertion) return

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted(i-1) <= DataUnsorted(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "LeaderUnsorted, LeaderSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted_RK(i), DataUnsorted(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        do i = 2, ndata
            assertion = assertion .and. DataUnsorted2(i-1) <= DataUnsorted2(i)
        end do

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))") 
            write(Test%outputUnit,"(*(g0,:,', '))") "RooterUnsorted, RooterSorted"
            do i = 1, ndata
                write(Test%outputUnit,"(*(g0,:,', '))") DataUnsorted2_RK(i), DataUnsorted2(i)
            end do
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_sortAscendingWithRooter_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMedian_RK_1() result(assertion)

        use Err_mod, only: Err_type
        implicit none

        logical                     :: assertion
        real(RK)    , parameter     :: median_ref = 5.477003450000000_RK
        real(RK)                    :: median
        type(Err_type)              :: Err

        call getMedian(lenArray = ndata, Array = DataUnsorted_RK, median = median, Err = Err)
        assertion = .not. Err%occurred
        if (.not. assertion) return

        assertion = median == median_ref

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))") 
            write(Test%outputUnit,"(*(g0,:,' '))") "median_ref  =", median_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "median      =", median
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMedian_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMedian_RK_2() result(assertion)

        use Err_mod, only: Err_type
        implicit none

        logical                     :: assertion
        real(RK)    , parameter     :: median_ref = 5.477415200000000_RK
        real(RK)                    :: median
        type(Err_type)              :: Err

        call getMedian_RK(lenArray = ndata-1, Array = DataUnsorted_RK(1:ndata-1), median = median, Err = Err)
        assertion = .not. Err%occurred
        if (.not. assertion) return

        assertion = median == median_ref

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))") 
            write(Test%outputUnit,"(*(g0,:,' '))") "median_ref  =", median_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "median      =", median
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMedian_RK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Sort_mod