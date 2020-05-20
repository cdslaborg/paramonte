!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module Sort_mod

    implicit none

    interface sortAscending
        module procedure :: sortAscending_RK, sortAscending2_RK
    end interface sortAscending

    interface indexArray
        module procedure :: indexArray_IK, indexArray_RK
    end interface indexArray

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    recursive subroutine sortArray(array)
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
  
!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine sortAscending_RK(np,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortAscending_RK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        implicit none
        integer(IK), parameter     :: nn = 15, nstack = 100
        integer(IK), intent(in)    :: np
        real(RK)   , intent(inout) :: Point(np)
        real(RK)                   :: dummy
        integer(IK)                :: k,i,j,jstack,m,r,istack(nstack)
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
                if (jstack > nstack) then
                    write(*,*) "sortAscending_RK() failed: nstack too small" ! xxx: needs improvement
                    error stop
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! Given Array(1:n) returns the array Indx(1:n) such that Array(Indx(j)), j=1:n is in ascending order.
    subroutine indexArray_RK(n,Array,Indx)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: indexArray_RK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        implicit none
        integer(IK), intent(in)  :: n
        real(RK)   , intent(in)  :: Array(n)
        integer(IK), intent(out) :: Indx(n)
        integer(IK), parameter   :: nn=15, nstack=50
        integer(IK)              :: k,i,j,indext,jstack,l,r
        integer(IK)              :: istack(nstack)
        real(RK)                 :: a
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
                if (jstack > nstack) then
                    write(*,*) "NSTACK too small in indexArray_RK()" ! xxx: needs improvement
                    error stop
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
        subroutine exchangeIndex(i,j)
            integer(IK), intent(inout) :: i,j
            integer(IK)                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
    end subroutine indexArray_RK

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! Given Array(1:n) returns the array Indx(1:n) such that Array(Indx(j)), j=1:n is in ascending order.
    subroutine indexArray_IK(n,Array,Indx)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: indexArray_IK
#endif
        use Constants_mod, only: IK, RK
        use Misc_mod, only: swap
        implicit none
        integer(IK), intent(in)  :: n
        integer(IK), intent(in)  :: Array(n)
        integer(IK), intent(out) :: Indx(n)
        integer(IK), parameter   :: nn=15_IK, nstack=50_IK
        integer(IK)              :: k,i,j,indext,jstack,l,r
        integer(IK)              :: istack(nstack)
        integer(IK)              :: a
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
                if (jstack > nstack) then
                    write(*,*) "NSTACK too small in indexArray_IK" ! xxx: needs improvement
                    stop
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
        subroutine exchangeIndex(i,j)
            integer(IK), intent(inout) :: i,j
            integer(IK)                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
  end subroutine indexArray_IK

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! Sorts Array(1:lenArray) in ascending order using Quicksort while making the corresponding rearrangement of Slave(1:lenArray).
    subroutine sortAscending2_RK(lenArray,Array,Slave)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortAscending2_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: lenArray
        real(RK), intent(inout) :: Array(lenArray), Slave(lenArray)
        integer(IK) :: Indx(lenArray)
        call indexArray_RK(lenArray,Array,Indx)
        Array = Array(Indx)
        Slave = Slave(Indx)
    end subroutine sortAscending2_RK

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Sort_mod