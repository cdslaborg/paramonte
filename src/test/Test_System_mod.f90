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

module Test_System_mod

    use System_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_System

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_System()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        if (Test%Image%isFirst) then
            call test_EnvVar_type()
            call test_CmdArg_type()
            call test_SystemInfo_type()
            call test_OS_type()
            call test_RandomFileName_type()
            call test_removeFile()
        end if
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_System

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_EnvVar_type()

        use String_mod, only: num2str
        implicit none
        type(EnvVar_type) :: EnvVar

        if (Test%Image%isFirst) call Test%testing("EnvVar_type")
        EnvVar%name = "OS"
        call EnvVar%get(EnvVar%name,EnvVar%value,EnvVar%length,EnvVar%Err)
        EnvVar%Err%msg = "Error occurred while querying OS type.\n" // EnvVar%Err%msg
        call Test%checkForErr(EnvVar%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "EnvVar%name     : ", EnvVar%name
            write(Test%outputUnit,"(2A)")   "EnvVar%value    : ", EnvVar%value
            write(Test%outputUnit,"(2A)")   "EnvVar%length   : ", num2str(EnvVar%length)
            write(Test%outputUnit,"(2A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_EnvVar_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_CmdArg_type()

        use String_mod, only: num2str
        implicit none
        integer           :: i
        type(CmdArg_type) :: CmdArg

        if (Test%Image%isFirst) call Test%testing("CmdArg_type")
        call CmdArg%query()
        CmdArg%Err%msg = "Error occurred while querying OS type.\n" // CmdArg%Err%msg
        call Test%checkForErr(CmdArg%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "CmdArg%cmd      : ", CmdArg%cmd
            write(Test%outputUnit,"(2A)")   "CmdArg%count    : ", num2str(CmdArg%count)
            write(Test%outputUnit,      "(*('CmdArg%slash    : ', 2A))") (CmdArg%Arg(i)%record, new_line('a'), i=1,CmdArg%count)
            write(Test%outputUnit,"(2A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_CmdArg_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_SystemInfo_type()

        implicit none
        integer           :: i
        type(SystemInfo_type) :: SystemInfo

        if (Test%Image%isFirst) call Test%testing("SystemInfo_type")

        SystemInfo = SystemInfo_type()
        call Test%checkForErr(SystemInfo%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(1A)")
            do i = 1, SystemInfo%nRecord
                write(Test%outputUnit,"(1A)") SystemInfo%List(i)%record
            end do
            write(Test%outputUnit,"(1A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_SystemInfo_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_OS_type()

        use String_mod, only: log2str
        implicit none
        type(OS_type) :: OS

        if (Test%Image%isFirst) call Test%testing("OS_type")
        call OS%query()
        OS%Err%msg = "Error occurred while querying OS type.\n" // OS%Err%msg
        call Test%checkForErr(OS%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "OS%name     : ", OS%name
            write(Test%outputUnit,"(2A)")   "OS%slash    : ", OS%slash
            write(Test%outputUnit,"(2A)")   "OS%isWindows: ", log2str(OS%isWindows)
            write(Test%outputUnit,"(2A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_OS_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_RandomFileName_type()

        use String_mod, only: log2str
        implicit none
        type(RandomFileName_type)   :: RFN
        logical                     :: fileExists

        if (Test%Image%isFirst) call Test%testing("RandomFileName_type")
        RFN = RandomFileName_type(key="test_RandomFileName_type")
        RFN%Err%msg = "Error occurred while generating random file name.\n" // RFN%Err%msg
        call Test%checkForErr(RFN%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "RFN%path     : ", RFN%path
            write(Test%outputUnit,"(2A)")   "RFN%dir      : ", RFN%dir
            write(Test%outputUnit,"(2A)")   "RFN%key      : ", RFN%key
            write(Test%outputUnit,"(2A)")   "RFN%ext      : ", RFN%ext
            write(Test%outputUnit,"(2A)")
        end if

        inquire(file=RFN%path,exist=fileExists)
        Test%assertion = .not.fileExists
        call Test%verify()
        Test%assertion = RFN%dir == ""
        call Test%verify()
        Test%assertion = RFN%key == "test_RandomFileName_type"  !  // "_image_" // num2str(this_image())
        call Test%verify()
        Test%assertion = RFN%ext == ".rfn"
        call Test%verify()

    end subroutine test_RandomFileName_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_removeFile()

        use String_mod, only: log2str
        implicit none
        type(RandomFileName_type)   :: RFN
        type(OS_type)               :: OS
        logical                     :: fileExists
        integer                     :: fileUnit

        if (Test%Image%isFirst) call Test%testing("removeFile()")
        RFN = RandomFileName_type(key="test_RandomFileName_type")
        RFN%Err%msg = "Error occurred while generating random file name.\n" // RFN%Err%msg
        call Test%checkForErr(RFN%Err)

        open(newunit=fileUnit,file=RFN%path,status="new")
        close(fileUnit)

        call OS%query()
        OS%Err%msg = "Error occurred while querying OS type.\n" // OS%Err%msg
        call Test%checkForErr(OS%Err)

        call removeFile(path=RFN%path,isWindows=OS%isWindows,Err=OS%Err)
        OS%Err%msg = "Error occurred.\n" // OS%Err%msg
        call Test%checkForErr(OS%Err)

        inquire(file=RFN%path,exist=fileExists)
        Test%assertion = .not.fileExists
        call Test%verify()

    end subroutine test_removeFile

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_System_mod
