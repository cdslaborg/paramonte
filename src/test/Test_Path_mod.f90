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

module Test_Path_mod

    use Path_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Path

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Path()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        
        if (Test%Image%isFirst) call test_Path_type()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_Path

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Path_type()

        use String_mod, only: num2str
        !use System_mod, only: executeCmd
        implicit none
        type(Path_type) :: Path


        call Test%testing("Path_type")

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        call Test%checkForErr(Path%Err)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%original  : ", Path%original
            write(Test%outputUnit,"(2A)")   "Path%modified  : ", Path%modified
            write(Test%outputUnit,"(2A)")   "Path%slashOS   : ", Path%slashOS
            write(Test%outputUnit,"(2A)")   "Path%dir       : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name      : ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext       : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%original == "./temp\ Folder/\{inside\}\/"
        call Test%verify()


        call Test%testing("Path_type%winify()")

        Path%original = "./temp\ Folder/\{inside\}/"
        call Path%winify(Path%original,Path%modified,Path%Err)
        if (Test%isDebugMode) then
            write(*,*) Path%original
            write(*,*) Path%modified
        end if
        call Test%checkForErr(Path%Err)
        Test%assertion = Path%modified == '".\temp Folder\{inside}\"'
        call Test%verify()


        call Test%testing("Path_type%linify()")

        Path%original = '".\temp Folder\{inside}\"'
        call Path%linify(Path%original,Path%modified)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%original: ", Path%original
            write(Test%outputUnit,"(2A)")   "Path%modified: ", Path%modified
            write(Test%outputUnit,"(2A)")   "Path%slashOS : ", Path%slashOS
            write(Test%outputUnit,"(2A)")   "Path%dir       : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name      : ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext       : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%modified == "./temp\ Folder/\{inside\}/"
        call Test%verify()


        call Test%testing("Path_type%linify()")

        Path%original = '".\temp Folder\{inside}\-"'
        call Path%linify(Path%original,Path%modified)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%original: ", Path%original
            write(Test%outputUnit,"(2A)")   "Path%modified: ", Path%modified
            write(Test%outputUnit,"(2A)")   "Path%slashOS : ", Path%slashOS
            write(Test%outputUnit,"(2A)")   "Path%dir       : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name      : ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext       : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%modified == "./temp\ Folder/\{inside\}/-"
        call Test%verify()


        call Test%testing("Path_type%getDirNameExt()")

        Path%original = ".\temp Folder\{inside}\"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%dir : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name: ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%dir == ".\temp Folder\{inside}\"
        call Test%verify()
        Test%assertion = Path%name == ""
        call Test%verify()
        Test%assertion = Path%ext == ""
        call Test%verify()


        call Test%testing("Path_type%getDirFullName()")

        Path%original = ".\temp Folder\{inside}\-"
        call Path%getDirFullName(Path%original,"\",Path%dir,Path%name)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%dir : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name: ", Path%name
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%dir == ".\temp Folder\{inside}\"
        call Test%verify()
        Test%assertion = Path%name == "-"
        call Test%verify()


        call Test%testing("Path_type%getDirFullName()")

        Path%original = "-"
        call Path%getNameExt(Path%original,Path%name,Path%ext)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%name: ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%name == "-"
        call Test%verify()
        Test%assertion = Path%ext == ""
        call Test%verify()

        Path%original = ".\temp Folder\{inside}\Temp.tXt"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "Path%dir : ", Path%dir
            write(Test%outputUnit,"(2A)")   "Path%name: ", Path%name
            write(Test%outputUnit,"(2A)")   "Path%ext : ", Path%ext
            write(Test%outputUnit,"(2A)")
        end if
        Test%assertion = Path%dir == ".\temp Folder\{inside}\"
        call Test%verify()
        Test%assertion = Path%name == "Temp"
        call Test%verify()
        Test%assertion = Path%ext == ".tXt"
        call Test%verify()

    end subroutine test_Path_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Path_mod
