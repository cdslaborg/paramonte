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

module Test_Misc_mod

    use Misc_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Misc

    type(Test_type) :: Test

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Misc()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_findUnique()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif
    end subroutine test_Misc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_findUnique()

        use Constants_mod, only: RK, IK

        implicit none
        integer(IK), parameter      :: VECTOR(*) = [1,2,1,3,5,5,2]
        integer(IK), parameter      :: LEN_VECTOR = size(VECTOR)
        integer(IK), parameter      :: UNIQUE_VALUE(*) = [1,2,3,5]
        integer(IK), parameter      :: UNIQUE_COUNT(*) = [2,2,1,2]
        integer(IK), allocatable    :: UniqueValue(:), UniqueCount(:), ZeroLenVector(:)
        integer(IK)                 :: lenUnique

        if (Test%Image%isFirst) call Test%testing("findUnique()")

        call findUnique ( lenVector = LEN_VECTOR &
                        , Vector = VECTOR &
                        , UniqueValue = UniqueValue &
                        , UniqueCount = UniqueCount &
                        , lenUnique = lenUnique &
                        )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", VECTOR
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_VALUE", UNIQUE_VALUE
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UNIQUE_COUNT", UNIQUE_COUNT
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = all(UniqueValue==UNIQUE_VALUE) .and. all(UniqueCount==UNIQUE_COUNT)
        call Test%verify()

        ! test with empty input vector

        allocate(ZeroLenVector(0))
        call findUnique ( lenVector = 0_IK &
                        , Vector = ZeroLenVector &
                        , UniqueValue = UniqueValue &
                        , UniqueCount = UniqueCount &
                        , lenUnique = lenUnique &
                        )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "VECTOR", ZeroLenVector
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueValue ", UniqueValue
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "UniqueCount ", UniqueCount
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "lenUnique", lenUnique
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        !Test%assertion = all(UniqueValue==UNIQUE_VALUE) .and. all(UniqueCount==UNIQUE_COUNT)
        !call Test%verify()
        !call Test%skipping()

    end subroutine test_findUnique

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Misc_mod