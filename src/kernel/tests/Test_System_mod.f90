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

!>  \brief This module contains tests of the module [System_mod](@ref system_mod).
!>  @author Amir Shahmoradi

module Test_System_mod

    use System_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_System

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_System()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        if (Test%Image%isFirst) then
            call Test%run(test_sleep_1, "test_sleep_1")
            call Test%run(test_OS_type_1, "test_OS_type_1")
            call Test%run(test_copyFile_1, "test_copyFile_1")
            call Test%run(test_removeFile_1, "test_removeFile_1")
            call Test%run(test_removeFile_2, "test_removeFile_2")
            call Test%run(test_executeCmd_1, "test_executeCmd_1")
            call Test%run(test_executeCmd_2, "test_executeCmd_2")
            call Test%run(test_EnvVar_type_1, "test_EnvVar_type_1")
            call Test%run(test_EnvVar_type_2, "test_EnvVar_type_2")
            call Test%run(test_EnvVar_type_3, "test_EnvVar_type_3")
            call Test%run(test_CmdArg_type_1, "test_CmdArg_type_1")
            call Test%run(test_SysCmd_type_1, "test_SysCmd_type_1")
            call Test%run(test_SystemInfo_type_1, "test_SystemInfo_type_1")
            call Test%run(test_RandomFileName_type_1, "test_RandomFileName_type_1")
        end if
        call Test%finalize()

    end subroutine test_System

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sleep_1() result(assertion)
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        type(Err_type)  :: Err
        call sleep(0.1_RK, Err)
        assertion = .not. Err%occurred
    end function test_sleep_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_SystemInfo_type_1() result(assertion)
        implicit none
        logical                 :: assertion
        type(SystemInfo_type)   :: SystemInfo
        SystemInfo = SystemInfo_type()
        assertion = .not. SystemInfo%Err%occurred .and. allocated(SystemInfo%List) .and. size(SystemInfo%List) > 0
    end function test_SystemInfo_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_EnvVar_type_1() result(assertion)

        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(EnvVar_type)       :: EnvVar

        EnvVar%name = ""
        call EnvVar%get(EnvVar%name,EnvVar%value,EnvVar%length,EnvVar%Err)
        assertion = EnvVar%Err%occurred

    end function test_EnvVar_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_EnvVar_type_2() result(assertion)

        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(EnvVar_type)       :: EnvVar

        EnvVar%name = "OS"
        call EnvVar%get(EnvVar%name,EnvVar%value,EnvVar%length,EnvVar%Err)
        assertion = .not. EnvVar%Err%occurred .and. allocated(EnvVar%name) .and. allocated(EnvVar%value)

        !if (Test%isDebugMode .and. .not. assertion) then
        !    write(Test%outputUnit,"(2A)")
        !    write(Test%outputUnit,"(2A)")   "EnvVar%name     : ", EnvVar%name
        !    write(Test%outputUnit,"(2A)")   "EnvVar%value    : ", EnvVar%value
        !    write(Test%outputUnit,"(2A)")   "EnvVar%length   : ", num2str(EnvVar%length)
        !    write(Test%outputUnit,"(2A)")
        !end if

    end function test_EnvVar_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_EnvVar_type_3() result(assertion)

        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(EnvVar_type)       :: EnvVar

        EnvVar%name = "OS"
        call EnvVar%get(EnvVar%name, EnvVar%value, Err = EnvVar%Err)
        assertion = .not. EnvVar%Err%occurred .and. allocated(EnvVar%name) .and. allocated(EnvVar%value)

    end function test_EnvVar_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_CmdArg_type_1() result(assertion)

        use String_mod, only: num2str
        implicit none
        !integer                 :: i
        logical                 :: assertion
        type(CmdArg_type)       :: CmdArg

        call CmdArg%query()
        assertion = .not. CmdArg%Err%occurred
        if (.not. assertion) return

        !if (Test%isDebugMode .and. .not. assertion) then
        !    write(Test%outputUnit,"(2A)")
        !    write(Test%outputUnit,"(2A)")   "CmdArg%cmd      : ", CmdArg%cmd
        !    write(Test%outputUnit,"(2A)")   "CmdArg%count    : ", num2str(CmdArg%count)
        !    write(Test%outputUnit,      "(*('CmdArg%slash    : ', 2A))") (CmdArg%Arg(i)%record, new_line('a'), i=1,CmdArg%count)
        !    write(Test%outputUnit,"(2A)")
        !end if

    end function test_CmdArg_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_OS_type_1() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                 :: assertion
        type(OS_type)           :: OS

        call OS%query(shellQueryEnabled = .true.)
        assertion = .not. OS%Err%occurred .and. .not. OS%Shell%Err%occurred
        if (.not. assertion) return

        !if (Test%isDebugMode .and. .not. assertion) then
        !    write(Test%outputUnit,"(2A)")
        !    write(Test%outputUnit,"(2A)")   "OS%name     : ", OS%name
        !    write(Test%outputUnit,"(2A)")   "OS%slash    : ", OS%slash
        !    write(Test%outputUnit,"(2A)")   "OS%isWindows: ", log2str(OS%isWindows)
        !    write(Test%outputUnit,"(2A)")
        !end if

    end function test_OS_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_SysCmd_type_1() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                 :: assertion
        type(SysCmd_type)       :: SysCmd
        type(OS_type)           :: OS

        call OS%query(shellQueryEnabled = .true.)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        if (OS%Shell%isUnix) then
            SysCmd = SysCmd_type("ls >/dev/null 2>&1", .true.)
            assertion = .not. SysCmd%Err%occurred
            if (.not. assertion) return
            SysCmd = SysCmd_type("ls >/dev/null 2>&1", .false.)
            assertion = .not. SysCmd%Err%occurred
            if (.not. assertion) return
        else
            SysCmd = SysCmd_type("dir > nul 2>&1", .true.)
            assertion = .not. SysCmd%Err%occurred
            if (.not. assertion) return
            SysCmd = SysCmd_type("dir > nul 2>&1", .false.)
            assertion = .not. SysCmd%Err%occurred
            if (.not. assertion) return
        end if

        !if (Test%isDebugMode .and. .not. assertion) then
        !    write(Test%outputUnit,"(2A)")
        !    write(Test%outputUnit,"(2A)")   "OS%name     : ", OS%name
        !    write(Test%outputUnit,"(2A)")   "OS%slash    : ", OS%slash
        !    write(Test%outputUnit,"(2A)")   "OS%isWindows: ", log2str(OS%isWindows)
        !    write(Test%outputUnit,"(2A)")
        !end if

    end function test_SysCmd_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_executeCmd_1() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(OS_type)               :: OS
        integer                     :: exitstat
        logical, parameter          :: wait = .true.
        character(:), allocatable   :: command

        call OS%query(shellQueryEnabled = .true.)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        if (OS%Shell%isUnix) then
            command = "ls >/dev/null 2>&1"
        else
            command = "dir > nul 2>&1"
        end if

        call executeCmd(command,wait,exitstat,OS%Err)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

    end function test_executeCmd_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_executeCmd_2() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(OS_type)               :: OS
        character(:), allocatable   :: command

        call OS%query(shellQueryEnabled = .true.)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        if (OS%Shell%isUnix) then
            command = "ls >/dev/null 2>&1"
        else
            command = "dir > nul 2>&1"
        end if

        call executeCmd(command)

    end function test_executeCmd_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_RandomFileName_type_1() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        logical                     :: fileExists

        RFN = RandomFileName_type(key="test_RandomFileName_type")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        inquire(file=RFN%path,exist=fileExists)
        assertion = assertion .and. .not. fileExists
        assertion = assertion .and. RFN%dir == ""
        assertion = assertion .and. RFN%key == "test_RandomFileName_type"  !  // "_image_" // num2str(this_image())
        assertion = assertion .and. RFN%ext == ".rfn"

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")   "RFN%path     : ", RFN%path
            write(Test%outputUnit,"(2A)")   "RFN%dir      : ", RFN%dir
            write(Test%outputUnit,"(2A)")   "RFN%key      : ", RFN%key
            write(Test%outputUnit,"(2A)")   "RFN%ext      : ", RFN%ext
            write(Test%outputUnit,"(2A)")
        end if

    end function test_RandomFileName_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_removeFile_1() result(assertion)

        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        logical                     :: fileExists
        integer                     :: fileUnit

        RFN = RandomFileName_type(key="test_RandomFileName_type")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        open(newunit=fileUnit,file=RFN%path,status="new")
        close(fileUnit)

        call removeFile(path=RFN%path,Err=RFN%Err)
        assertion = assertion .and. .not. RFN%Err%occurred
        if (.not. assertion) return

        inquire(file=RFN%path,exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_removeFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_removeFile_2() result(assertion)

        use Constants_mod, only: RK
        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        logical                     :: fileExists
        integer                     :: fileUnit

        RFN = RandomFileName_type(key="test_RandomFileName_type")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        open(newunit=fileUnit,file=RFN%path,status="new")
        close(fileUnit)

        call removeFile(path=RFN%path)
        inquire(file=RFN%path,exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_removeFile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyFile_1() result(assertion)

        use Constants_mod, only: RK
        use String_mod, only: log2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        type(OS_type)               :: OS
        logical                     :: fileExists
        integer                     :: fileUnit

        RFN = RandomFileName_type(key="test_RandomFileName_type")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        open(newunit=fileUnit,file=RFN%path,status="new")
        close(fileUnit)

        call OS%query(shellQueryEnabled = .true.)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        call copyFile(RFN%path, RFN%path//".copy", .not. OS%Shell%isUnix, OS%Err)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        inquire(file=RFN%path//".copy",exist=fileExists)
        assertion = assertion .and. fileExists

        call removeFile(RFN%path)
        inquire(file=RFN%path,exist=fileExists)
        assertion = assertion .and. .not. fileExists

        call removeFile(RFN%path//".copy")
        inquire(file=RFN%path//".copy",exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_copyFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_System_mod
