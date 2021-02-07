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

!>  \brief This module contains miscellaneous procedures.
!>  \author Amir Shahmoradi

module Unique_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Unique_mod"

    interface findUnique
        module procedure :: findUnique_IK
    end interface findUnique

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Find the unique values in the input integer vector.
    !>
    !> @param[in]       lenVector       :   The size of the input square matrix - `nd` by `nd`.
    !> @param[in]       Vector          :   The input integer vector.
    !> @param[out]      lenUnique       :   The length of `UniqueValue`, that is, the total number of unique values.
    !> @param[out]      UniqueValue     :   The vector of unique values identified in the input vector.
    !> @param[out]      UniqueCount     :   The counts of each unique value in the input vector.
    !> @param[out]      UniqueIndex     :   A jaggedArray of type [IntVec_type](@ref jaggedarray_mod::intvec_type) of length `lenUnique`,
    !>                                      the ith element of which is a vector of length `UniqueCount(i)` that contains the
    !>                                      indices of `Vector` where `UniqueValue(i)` occur (**optional**).
    !> @param[out]      Err             :   An object of type [Err_type](@ref err_mod::err_type). If present, the output `UniqueCount` 
    !>                                      will be sorted **descending** and along with it `UniqueValue` and `UniqueIndex` (**optional**).
    !> \warning
    !> To avoid extra data copy and improve performance, the output arrays `UniqueCount` and `UniqueValue` will not be
    !> resized from `(1:lenVector)` to `(1:lenUnique)`. The onus is on the user to ensure only the elements `(1:lenUnique)`
    !> are used for any subsequent work as only these elements are meaningful. **However**, if `Err` argument is present,
    !> ther two aforementioned arrays will be ordered and then automatically resized to `(1:lenUnique)`.
    pure subroutine findUnique_IK(lenVector, Vector, lenUnique, UniqueValue, UniqueCount, UniqueIndex, Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: findUnique_IK
#endif
        use JaggedArray_mod, only: IV => IntVec_type
        use Sort_mod, only: indexArray
        use Constants_mod, only: IK
        use Err_mod, only: Err_type
        implicit none
        integer(IK)     , intent(in)                                :: lenVector
        integer(IK)     , intent(in)                                :: Vector(lenVector)
        integer(IK)     , intent(out)                               :: lenUnique
        integer(IK)     , intent(out)   , allocatable               :: UniqueValue(:)
        integer(IK)     , intent(out)   , allocatable               :: UniqueCount(:)
        type(IV)        , intent(out)   , allocatable   , optional  :: UniqueIndex(:)
        type(Err_type)  , intent(out)                   , optional  :: Err
        integer(IK)                     , allocatable               :: Indx(:)
        integer(IK)                                                 :: ivec, iuniq, counter
        logical                                                     :: isUnique
        allocate(UniqueValue(lenVector))
        allocate(UniqueCount(lenVector), source = 0_IK)

        lenUnique = 0

        do ivec = 1, lenVector
            isUnique = .true.
            loopSearchUnique: do iuniq = 1, lenUnique
                if (UniqueValue(iuniq)==Vector(ivec)) then
                    UniqueCount(iuniq) = UniqueCount(iuniq) + 1
                    isUnique = .false.
                    exit loopSearchUnique
                end if
            end do loopSearchUnique
            if (isUnique) then
                lenUnique = lenUnique + 1
                UniqueValue(lenUnique) = Vector(ivec)
                UniqueCount(lenUnique) = UniqueCount(lenUnique) + 1
            end if
        end do

        if (present(Err)) then
            allocate(Indx(lenUnique))
            call indexArray(lenUnique, UniqueCount(1:lenUnique), Indx, Err)
            if (Err%occurred) return
            UniqueCount = UniqueCount(Indx(lenUnique:1:-1))
            UniqueValue = UniqueValue(Indx(lenUnique:1:-1))
        end if

        if (present(UniqueIndex)) then
            allocate(UniqueIndex(lenUnique))
            loopUniqueIndex: do iuniq = 1, lenUnique
                allocate(UniqueIndex(iuniq)%Vector(UniqueCount(iuniq)))
                counter = 1_IK
                loopOverVector: do ivec = 1, lenVector
                    if (UniqueValue(iuniq) == Vector(ivec)) then
                        UniqueIndex(iuniq)%Vector(counter) = ivec
                        counter = counter + 1_IK
                    end if
                    if (counter > UniqueCount(iuniq)) exit loopOverVector
                end do loopOverVector
            end do loopUniqueIndex
        end if

    end subroutine findUnique_IK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Unique_mod