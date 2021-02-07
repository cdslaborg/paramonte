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

!>  \brief This module contains tests of the module [Unique_mod](@ref misc_mod).
!>  \author Amir Shahmoradi

module Test_Unique_mod

    use Unique_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Unique

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Unique()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_findUnique_1, "test_findUnique_1")
        call Test%finalize()
    end subroutine test_Unique

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_findUnique_1() result(assertion)

        use Constants_mod, only: RK, IK
        use JaggedArray_mod, only: IV => IntVec_type

        implicit none
        logical                     :: assertion
        logical                     :: assertionCurrent
        integer(IK) , parameter     :: VECTOR(*) = [1,2,1,3,5,5,2]
        integer(IK) , parameter     :: LEN_VECTOR = size(VECTOR)
        integer(IK) , parameter     :: UNIQUE_VALUE(*) = [1,2,3,5]
        integer(IK) , parameter     :: UNIQUE_COUNT(*) = [2,2,1,2]
        integer(IK) , allocatable   :: ZeroLenVector(:)
        integer(IK) , allocatable   :: UniqueValue(:)
        integer(IK) , allocatable   :: UniqueCount(:)
        type(IV)    , allocatable   :: UniqueIndex(:)
        type(IV)    , allocatable   :: UNIQUE_INDEX(:)
        integer(IK)                 :: lenUnique, i
        type(Err_type)              :: Err

        call findUnique ( lenVector = LEN_VECTOR &
                        , Vector = VECTOR &
                        , lenUnique = lenUnique &
                        , UniqueValue = UniqueValue &
                        , UniqueCount = UniqueCount &
                        )

        assertion = all(UniqueValue(1:lenUnique)==UNIQUE_VALUE) .and. all(UniqueCount(1:lenUnique)==UNIQUE_COUNT)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", VECTOR
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! test UniqueIndex

        allocate(UNIQUE_INDEX(size(UNIQUE_COUNT)))
        UNIQUE_INDEX(1)%Vector = [5,6]
        UNIQUE_INDEX(2)%Vector = [2,7]
        UNIQUE_INDEX(3)%Vector = [1,3]
        UNIQUE_INDEX(4)%Vector = [4]

        call findUnique ( lenVector = LEN_VECTOR &
                        , Vector = VECTOR &
                        , lenUnique = lenUnique &
                        , UniqueValue = UniqueValue &
                        , UniqueCount = UniqueCount &
                        , UniqueIndex = UniqueIndex &
                        , Err = Err &
                        )

        assertion = assertion .and. .not. Err%occurred

        do i = 1, lenUnique

            assertionCurrent = all(UniqueIndex(i)%Vector == UNIQUE_INDEX(i)%Vector)
            assertion = assertion .and. assertionCurrent

            if (Test%isVerboseMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
                write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP

            if (i>1_IK) assertionCurrent = assertionCurrent .and. UniqueCount(i) <= UniqueCount(i-1)
            assertion = assertion .and. assertionCurrent

            if (Test%isVerboseMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", VECTOR
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
                write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
                write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_INDEX(i)%Vector ", UNIQUE_INDEX(i)%Vector
                write(Test%outputUnit,"(*(g0,:,', '))") "UniqueIndex(i)%Vector  ", UniqueIndex(i)%Vector
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            ! LCOV_EXCL_STOP

            if (.not. assertion) return ! LCOV_EXCL_LINE

        end do

        ! test with empty input vector

        allocate(ZeroLenVector(0))
        call findUnique ( lenVector = 0_IK & ! LCOV_EXCL_LINE
                        , Vector = ZeroLenVector & ! LCOV_EXCL_LINE
                        , UniqueValue = UniqueValue & ! LCOV_EXCL_LINE
                        , UniqueCount = UniqueCount & ! LCOV_EXCL_LINE
                        , lenUnique = lenUnique & ! LCOV_EXCL_LINE
                        )

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", VECTOR
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", ZeroLenVector
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount(1:lenUnique)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_findUnique_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Unique_mod ! LCOV_EXCL_LINE