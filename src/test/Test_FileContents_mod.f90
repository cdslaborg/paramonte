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
