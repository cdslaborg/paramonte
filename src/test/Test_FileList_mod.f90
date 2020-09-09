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

module Test_FileList_mod

    use FileList_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_FileList

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_FileList()

        Test = Test_type(moduleName=MODULE_NAME)
        if (Test%Image%isFirst) call test_FileList_type()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif

        
    end subroutine test_FileList

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_FileList_type()
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        type(FileList_type) :: FileList
        integer(IK)         :: i

        call Test%testing("FileList_type")

        if (Test%Image%isFirst) then
            FileList = FileList_type()
            call Test%checkForErr(FileList%Err)
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "searchStr  : ", FileList%searchStr
            write(Test%outputUnit,"(2A)")   "orderStr   : ", FileList%orderStr
            write(Test%outputUnit,"(2A)")   "excludeStr : ", FileList%excludeStr
            do i = 1,FileList%count
                write(Test%outputUnit,"(2A)")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(2A)")
        end if

        if (Test%Image%isFirst) then
            FileList = FileList_type(orderStr="name")
            call Test%checkForErr(FileList%Err)
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "searchStr  : ", FileList%searchStr
            write(Test%outputUnit,"(2A)")   "orderStr   : ", FileList%orderStr
            write(Test%outputUnit,"(2A)")   "excludeStr : ", FileList%excludeStr
            do i = 1,FileList%count
                write(Test%outputUnit,"(2A)")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(2A)")
        end if

        if (Test%Image%isFirst) then
            FileList = FileList_type(orderStr="date")
            call Test%checkForErr(FileList%Err)
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "searchStr  : ", FileList%searchStr
            write(Test%outputUnit,"(2A)")   "orderStr   : ", FileList%orderStr
            write(Test%outputUnit,"(2A)")   "excludeStr : ", FileList%excludeStr
            do i = 1,FileList%count
                write(Test%outputUnit,"(2A)")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(2A)")
        end if

        if (Test%Image%isFirst) then
            FileList = FileList_type(orderStr="name")
            call Test%checkForErr(FileList%Err)
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "searchStr  : ", FileList%searchStr
            write(Test%outputUnit,"(2A)")   "orderStr   : ", FileList%orderStr
            write(Test%outputUnit,"(2A)")   "excludeStr : ", FileList%excludeStr
            do i = 1,FileList%count
                write(Test%outputUnit,"(2A)")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(2A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_FileList_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_FileList_mod
