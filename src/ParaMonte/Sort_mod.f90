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

!>  \brief This module contains procedures for sorting arrays.
!>  @author Amir Shahmoradi

module Sort_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Sort_mod"

    interface sortAscending
        module procedure :: sortAscending_RK, sortAscending2_RK
    end interface sortAscending

    interface indexArray
        module procedure :: indexArray_IK, indexArray_RK
    end interface indexArray

    interface median
        module procedure :: median_RK
    end interface median

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Sort (recursively) the input real array of arbitrary size from the smallest value to the largest and
    !> return the result inside the input array.
    !>
    !> @param[inout] array : The input vector to be sorted.
    !>
    !> \warning
    !> On return, the contents of the input array is completely overwritten.
    pure recursive subroutine sortArray(array)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortArray
#endif
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(inout) :: array(:)
        integer(IK)             :: iq

        if(size(array) > 1) then
            call partition(array, iq)
            call sortArray(array(:iq-1))
            call sortArray(array(iq:))
        endif
    end subroutine sortArray

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine partition(array, marker)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: partition
#endif
        use Constants_mod, only: IK, RK
        implicit none
        real(RK)   , intent(inout) :: array(:)
        integer(IK), intent(out)   :: marker
        integer(IK)                :: i, j
        real(RK)                   :: temp
        real(RK)                   :: x      ! pivot point
        x = array(1)
        i= 0
        j= size(array) + 1
        do
            j = j-1
            do
                if (array(j) <= x) exit
                j = j-1
            end do
            i = i+1
            do
                if (array(i) >= x) exit
                i = i+1
            end do
            if (i < j) then ! exchange array(i) and array(j)
                temp = array(i)
                array(i) = array(j)
                array(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do
  end subroutine partition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Sort (recursively) the input real array of arbitrary size from the smallest value to the largest and
    !> return the result inside the input array.
    !>
    !> @param[in]       np      :   The length of the input vector to be sorted.
    !> @param[inout]    Point   :   The input vector to be sorted.
    !> @param[out]      Err     :   An object of class [Err_type](@ref err_mod::err_type).
    !>
    !> \warning
    !> On return, the contents of the input array is completely overwritten by the output sorted array.
    !>
    !> \warning
    !> On return, the value of `Err%%occurred` must be checked for any potential occurrences of errors during sorting.
    pure subroutine sortAscending_RK(np,Point,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortAscending_RK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , parameter     :: nn = 15, NSTACK = 100
        integer(IK)     , intent(in)    :: np
        real(RK)        , intent(inout) :: Point(np)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@sortAscending_RK()"
        real(RK)                        :: dummy
        integer(IK)                     :: k,i,j,jstack,m,r,istack(NSTACK)

        Err%occurred = .false.

        jstack=0
        m = 1
        r = np
        do
            if (r-m < nn) then
                do j = m+1,r
                    dummy = Point(j)
                    do i = j-1,m,-1
                        if (Point(i) <=  dummy) exit
                        Point(i+1) = Point(i)
                    end do
                    Point(i+1) = dummy
                end do
                if (jstack == 0) return
                r = istack(jstack)
                m = istack(jstack-1)
                jstack = jstack-2
            else
                k = (m+r)/2
                call swap(Point(k),Point(m+1))
                if (Point(m)>Point(r)) call swap(Point(m),Point(r))
                if (Point(m+1)>Point(r)) call swap(Point(m+1),Point(r))
                if (Point(m)>Point(m+1)) call swap(Point(m),Point(m+1))
                i = m+1
                j = r
                dummy = Point(m+1)
                do
                    do
                        i = i+1
                        if (Point(i) >= dummy) exit
                    end do
                    do
                        j = j-1
                        if (Point(j) <= dummy) exit
                    end do
                    if (j < i) exit
                    call swap(Point(i),Point(j))
                end do
                Point(m+1) = Point(j)
                Point(j) = dummy
                jstack = jstack+2
                if (jstack > NSTACK) then
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": NSTACK is too small."
                    return
                end if
                if (r-i+1 >= j-m) then
                    istack(jstack) = r
                    istack(jstack-1) = i
                    r = j-1
                else
                    istack(jstack) = j-1
                    istack(jstack-1) = m
                    m = i
                end if
            end if
        end do
    end subroutine sortAscending_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given the real `Array(1:n)` return the array `Indx(1:n)` such that `Array(Indx(j)), j=1:n` is in ascending order.
    !>
    !> @param[in]   n       :   The length of the input vector to be sorted.
    !> @param[in]   Array   :   The input vector of length `n` to be sorted.
    !> @param[out]  Indx    :   The output integer vector of indices of the sorted vector.
    !> @param[out]  Err     :   An object of class [Err_type](@ref err_mod::err_type).
    !>
    !> \warning
    !> On return, the value of `Err%%occurred` must be checked for any potential occurrences of errors during sorting.
    pure subroutine indexArray_RK(n,Array,Indx,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: indexArray_RK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: n
        real(RK)        , intent(in)    :: Array(n)
        integer(IK)     , intent(out)   :: Indx(n)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@indexArray_RK()"
        integer(IK)     , parameter     :: nn=15, NSTACK=50
        integer(IK)                     :: k,i,j,indext,jstack,l,r
        integer(IK)                     :: istack(NSTACK)
        real(RK)                        :: a

        Err%occurred = .false.

        do j = 1,n
            Indx(j) = j
        end do
        jstack=0
        l=1
        r=n
        do
            if (r-l < nn) then
                do j=l+1,r
                    indext=Indx(j)
                    a=Array(indext)
                    do i=j-1,l,-1
                        if (Array(Indx(i)) <= a) exit
                        Indx(i+1)=Indx(i)
                    end do
                    Indx(i+1)=indext
                end do
                if (jstack == 0) return
                r=istack(jstack)
                l=istack(jstack-1)
                jstack=jstack-2
            else
                k=(l+r)/2
                call swap(Indx(k),Indx(l+1))
                call exchangeIndex(Indx(l),Indx(r))
                call exchangeIndex(Indx(l+1),Indx(r))
                call exchangeIndex(Indx(l),Indx(l+1))
                i=l+1
                j=r
                indext=Indx(l+1)
                a=Array(indext)
                do
                    do
                        i=i+1
                        if (Array(Indx(i)) >= a) exit
                    end do
                    do
                        j=j-1
                        if (Array(Indx(j)) <= a) exit
                    end do
                    if (j < i) exit
                    call swap(Indx(i),Indx(j))
                end do
                Indx(l+1)=Indx(j)
                Indx(j)=indext
                jstack=jstack+2
                if (jstack > NSTACK) then
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": NSTACK is too small."
                    return
                end if
                if (r-i+1 >= j-l) then
                    istack(jstack)=r
                    istack(jstack-1)=i
                    r=j-1
                else
                    istack(jstack)=j-1
                    istack(jstack-1)=l
                    l=i
                end if
            end if
        end do
    contains
        pure subroutine exchangeIndex(i,j)
            integer(IK), intent(inout) :: i,j
            integer(IK)                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
    end subroutine indexArray_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given the integer `Array(1:n)` return the array `Indx(1:n)` such that `Array(Indx(j)), j=1:n` is in ascending order.
    !>
    !> @param[in]   n       :   The length of the input vector to be sorted.
    !> @param[in]   Array   :   The input vector of length `n` to be sorted.
    !> @param[out]  Indx    :   The output integer vector of indices of the sorted vector.
    !> @param[out]  Err     :   An object of class [Err_type](@ref err_mod::err_type).
    !>
    !> \warning
    !> On return, the value of `Err%%occurred` must be checked for any potential occurrences of errors during sorting.
    pure subroutine indexArray_IK(n,Array,Indx,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: indexArray_IK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: n
        integer(IK)     , intent(in)    :: Array(n)
        integer(IK)     , intent(out)   :: Indx(n)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@indexArray_IK()"
        integer(IK)     , parameter     :: nn=15_IK, NSTACK=50_IK
        integer(IK)                     :: k,i,j,indext,jstack,l,r
        integer(IK)                     :: istack(NSTACK)
        integer(IK)                     :: a

        Err%occurred = .false.

        do j = 1,n
            Indx(j) = j
        end do
        jstack=0
        l=1
        r=n
        do
            if (r-l < nn) then
                do j=l+1,r
                    indext=Indx(j)
                    a=Array(indext)
                    do i=j-1,l,-1
                        if (Array(Indx(i)) <= a) exit
                        Indx(i+1)=Indx(i)
                    end do
                    Indx(i+1)=indext
                end do
                if (jstack == 0) return
                r=istack(jstack)
                l=istack(jstack-1)
                jstack=jstack-2
            else
                k=(l+r)/2
                call swap(Indx(k),Indx(l+1))
                call exchangeIndex(Indx(l),Indx(r))
                call exchangeIndex(Indx(l+1),Indx(r))
                call exchangeIndex(Indx(l),Indx(l+1))
                i=l+1
                j=r
                indext=Indx(l+1)
                a=Array(indext)
                do
                    do
                        i=i+1
                        if (Array(Indx(i)) >= a) exit
                    end do
                    do
                        j=j-1
                        if (Array(Indx(j)) <= a) exit
                    end do
                    if (j < i) exit
                    call swap(Indx(i),Indx(j))
                end do
                Indx(l+1)=Indx(j)
                Indx(j)=indext
                jstack=jstack+2
                if (jstack > NSTACK) then
                    Err%occurred = .true.
                    Err%msg = PROCEDURE_NAME//": NSTACK is too small."
                    return
                end if
                if (r-i+1 >= j-l) then
                    istack(jstack)=r
                    istack(jstack-1)=i
                    r=j-1
                else
                    istack(jstack)=j-1
                    istack(jstack-1)=l
                    l=i
                end if
            end if
        end do
    contains
        pure subroutine exchangeIndex(i,j)
            integer(IK), intent(inout) :: i,j
            integer(IK)                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
  end subroutine indexArray_IK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Sort the real `Array(1:lenArray)` in ascending order using Quicksort while making the corresponding rearrangement of real `Slave(1:lenArray)`.
    !>
    !> @param[in]       lenArray    :   The length of the input vector to be sorted.
    !> @param[inout]    Array       :   The vector of length `lenArray` to be sorted.
    !> @param[inout]    Slave       :   The vector of length `lenArray` to be sorted according to the rearrangement of the elements of `Array`.
    !> @param[out]      Err         :   An object of class [Err_type](@ref err_mod::err_type).
    !>
    !> \warning
    !> On return, the value of `Err%%occurred` must be checked for any potential occurrences of errors during sorting.
    pure subroutine sortAscending2_RK(lenArray,Array,Slave,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortAscending2_RK
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: lenArray
        real(RK)        , intent(inout) :: Array(lenArray), Slave(lenArray)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@indexArray_IK()"
        integer(IK)                     :: Indx(lenArray)

        call indexArray_RK(lenArray,Array,Indx,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME//": NSTACK is too small."
            return
        end if
        Array = Array(Indx)
        Slave = Slave(Indx)
    end subroutine sortAscending2_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the median of the input vector.
    !>
    !> @param[in]       lenArray    :   The length of the input vector.
    !> @param[in]       Array       :   The input vector of length `lenArray` whose median is to be found.
    !> @param[out]      median      :   The median of the input `Array`.
    !> @param[out]      Err         :   An object of class [Err_type](@ref err_mod::err_type).
    !>
    !> \warning
    !> On return, the value of `Err%%occurred` must be checked for any potential occurrences of errors during sorting.
    pure subroutine median_RK(lenArray,Array,median,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: median_RK
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: lenArray
        real(RK)        , intent(in)    :: Array(lenArray)
        real(RK)        , intent(out)   :: median
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@median_RK()"
        real(RK)                        :: ArrayDummy(lenArray)

        ArrayDummy = Array
        call sortAscending(np=lenArray,Point=ArrayDummy,Err=Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME//Err%msg
            return
        end if
        median = ArrayDummy(lenArray/2_IK)

    end subroutine median_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Sort_mod