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
