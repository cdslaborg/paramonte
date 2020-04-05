!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
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

module Test_FileContents_mod

    use FileContents_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_FileContents

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_FileContents()

        Test = Test_type(moduleName=MODULE_NAME)
        if (Test%Image%isFirst) call test_FileContents_type()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_FileContents

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_FileContents_type()
        use System_mod, only: RandomFileName_type, removeFile, OS_type
        use String_mod, only: num2str
        implicit none
        type(RandomFileName_type)   :: RandomFileName
        type(FileContents_type)     :: FileContents
        type(OS_type)               :: OS
        integer                     :: fileUnit,i

        call Test%testing("FileContents_type")

        RandomFileName = RandomFileName_type(key="test_FileContents_type")
        call Test%checkForErr(RandomFileName%Err)

        open(newunit=fileUnit,file=RandomFileName%path,status="new")
        write(fileUnit,"(A)") "Testing FileContents_type..."
        write(fileUnit,"(A)") " "
        write(fileUnit,"(A)") ""
        write(fileUnit,"(A)")
        FileContents = FileContents_type(RandomFileName%path)
        call Test%checkForErr(FileContents%Err)
        call OS%query()
        call Test%checkForErr(OS%Err)
        call removeFile(RandomFileName%path,OS%isWindows,OS%Err)
        call Test%checkForErr(OS%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "numRecord  : ", num2str(FileContents%numRecord)
            do i = 1,FileContents%numRecord
                write(Test%outputUnit,"(2A)")   "FileContents%Line(" // num2str(i) // ")%record : ", FileContents%Line(i)%record
            end do
            write(Test%outputUnit,"(2A)")
        end if

        Test%assertion = FileContents%numRecord==4
        call Test%verify()

        Test%assertion = FileContents%Line(1)%record=="Testing FileContents_type..."
        call Test%verify()
        Test%assertion = FileContents%Line(2)%record==""
        call Test%verify()
        Test%assertion = FileContents%Line(3)%record==""
        call Test%verify()
        Test%assertion = FileContents%Line(4)%record==""
        call Test%verify()

    end subroutine test_FileContents_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_FileContents_mod
