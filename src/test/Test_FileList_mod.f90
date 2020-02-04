!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

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
