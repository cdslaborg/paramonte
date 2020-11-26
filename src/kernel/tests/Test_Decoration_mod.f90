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

!>  \brief This module contains tests of the module [Decoration_mod](@ref decoration_mod).
!>  @author Amir Shahmoradi

module Test_Decoration_mod

    !use, intrinsic :: iso_fortran_env, only: output_unit
    use Test_mod, only: Test_type
    use Constants_mod, only: IK
    use Decoration_mod

    implicit none

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Decoration()

        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_wrapText, "test_wrapText")
        call Test%run(test_drawLine_1, "test_drawLine_1")
        call Test%run(test_drawLine_2, "test_drawLine_2")
        call Test%run(test_drawLine_3, "test_drawLine_3")
        call Test%run(test_sandwich_1, "test_sandwich_1")
        call Test%run(test_sandwich_2, "test_sandwich_2")
        call Test%run(test_sandwich_3, "test_sandwich_3")
        call Test%run(test_sandwich_4, "test_sandwich_4")
        call Test%run(test_sandwich_5, "test_sandwich_5")
        call Test%run(test_getGenericFormat_1, "test_getGenericFormat_1")
        call Test%run(test_getGenericFormat_2, "test_getGenericFormat_2")
        call Test%run(test_getGenericFormat_3, "test_getGenericFormat_3")
        call Test%run(test_getGenericFormat_4, "test_getGenericFormat_4")
        call Test%run(test_getGenericFormat_5, "test_getGenericFormat_5")
        call Test%run(test_writeDecoratedText_1, "test_writeDecoratedText_1")
        call Test%run(test_writeDecoratedText_2, "test_writeDecoratedText_2")
        call Test%run(test_writeDecoratedList_1, "test_writeDecoratedList_1")
        call Test%run(test_writeDecoratedList_2, "test_writeDecoratedList_2")
        call Test%finalize()

    end subroutine test_Decoration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getGenericFormat_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: genericFormat_ref = "('ParaMonte',*(g25.10,:,','))"
        character(:), allocatable   :: genericFormat
        genericFormat = getGenericFormat(width = 25_IK, precision = 10_IK, delim = ",", prefix = "ParaMonte")
        assertion = genericFormat == genericFormat_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat_ref    =", genericFormat_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat        =", genericFormat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_getGenericFormat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getGenericFormat_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: genericFormat_ref = "(*(g25.10,:,','))"
        character(:), allocatable   :: genericFormat
        genericFormat = getGenericFormat(width = 25_IK, precision = 10_IK, delim = ",")
        assertion = genericFormat == genericFormat_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat_ref    =", genericFormat_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat        =", genericFormat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_getGenericFormat_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getGenericFormat_3() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: genericFormat_ref = "(*(g25.10))"
        character(:), allocatable   :: genericFormat
        genericFormat = getGenericFormat(width = 25_IK, precision = 10_IK)
        assertion = genericFormat == genericFormat_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat_ref    =", genericFormat_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat        =", genericFormat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_getGenericFormat_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getGenericFormat_4() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: genericFormat_ref = "(*(g25))"
        character(:), allocatable   :: genericFormat
        genericFormat = getGenericFormat(width = 25_IK)
        assertion = genericFormat == genericFormat_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat_ref    =", genericFormat_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat        =", genericFormat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_getGenericFormat_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getGenericFormat_5() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: genericFormat_ref = "(*(g0))"
        character(:), allocatable   :: genericFormat
        genericFormat = getGenericFormat()
        assertion = genericFormat == genericFormat_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat_ref    =", genericFormat_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "genericFormat        =", genericFormat
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_getGenericFormat_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_writeDecoratedText_1() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical                         :: assertion
        logical                         :: assertionCurrent
        type(Decoration_type)           :: Decoration
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        integer(IK)                     :: fileUnit, i, iostat
        integer(IK), parameter          :: NLINE = 19_IK

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref( 1)%record = ""
        OutputList_ref( 2)%record = ""
        OutputList_ref( 3)%record = "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        OutputList_ref( 4)%record = "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        OutputList_ref( 5)%record = "&&&&                        &&&&"
        OutputList_ref( 6)%record = "&&&&                        &&&&"
        OutputList_ref( 7)%record = "&&&&Have you asked yourself:&&&&"
        OutputList_ref( 8)%record = "&&&&                        &&&&"
        OutputList_ref( 9)%record = "&&&&s the Universe bother to&&&&"
        OutputList_ref(10)%record = "&&&&                        &&&&"
        OutputList_ref(11)%record = "&&&& the origin of mass and &&&&"
        OutputList_ref(12)%record = "&&&&                        &&&&"
        OutputList_ref(13)%record = "&&&&at is the origin of life&&&&"
        OutputList_ref(14)%record = "&&&&                        &&&&"
        OutputList_ref(15)%record = "&&&&                        &&&&"
        OutputList_ref(16)%record = "&&&&                        &&&&"
        OutputList_ref(17)%record = "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        OutputList_ref(18)%record = "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        OutputList_ref(19)%record = ""

        allocate(Decoration%List(4))

        Decoration%List(1)%record = "Have you asked yourself:"
        Decoration%List(2)%record = "Why does the Universe bother to exist?"
        Decoration%List(3)%record = "What is the origin of mass and matter?"
        Decoration%List(4)%record = "What is the origin of life?"

        Decoration%text =   "\n\n" // &
                            Decoration%List(1)%record // "\n\n" // &
                            Decoration%List(2)%record // "\n\n" // &
                            Decoration%List(3)%record // "\n\n" // &
                            Decoration%List(4)%record // "\n\n\n"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"Test_Decoration_mod@test_writeDecoratedText_1."//num2str(Test%Image%id)//".out", status = "replace")

        call Decoration%writeDecoratedText  ( Decoration%text &
                                            , newLine="\n" &
                                            , width = 32_IK &
                                            , symbol = "&" &
                                            , thicknessHorz = 4_IK &
                                            , thicknessVert = 2_IK &
                                            , marginTop = 2_IK &
                                            , marginBot = 1_IK &
                                            , outputUnit = fileUnit &
                                            )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)
            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            if (iostat/=0_IK) then
                assertion = .false.
                return
            end if
            OutputList(i)%record = trim(adjustl(OutputList(i)%record))

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if

        end do

    end function test_writeDecoratedText_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_writeDecoratedText_2() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK, NLC
        use String_mod, only: num2str
        implicit none
        logical                         :: assertion
        logical                         :: assertionCurrent
        type(Decoration_type)           :: Decoration
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        integer(IK)                     :: fileUnit, i, iostat
        integer(IK), parameter          :: NLINE = 14_IK

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref( 1)%record = "************************************************************************************************************************************"
        OutputList_ref( 2)%record = "****                                                                                                                            ****"
        OutputList_ref( 3)%record = "****                                                                                                                            ****"
        OutputList_ref( 4)%record = "****                                                  Have you asked yourself:                                                  ****"
        OutputList_ref( 5)%record = "****                                                                                                                            ****"
        OutputList_ref( 6)%record = "****                                           Why does the Universe bother to exist?                                           ****"
        OutputList_ref( 7)%record = "****                                                                                                                            ****"
        OutputList_ref( 8)%record = "****                                           What is the origin of mass and matter?                                           ****"
        OutputList_ref( 9)%record = "****                                                                                                                            ****"
        OutputList_ref(10)%record = "****                                                What is the origin of life?                                                 ****"
        OutputList_ref(11)%record = "****                                                                                                                            ****"
        OutputList_ref(12)%record = "****                                                                                                                            ****"
        OutputList_ref(13)%record = "****                                                                                                                            ****"
        OutputList_ref(14)%record = "************************************************************************************************************************************"

        allocate(Decoration%List(4))

        Decoration%List(1)%record = "Have you asked yourself:"
        Decoration%List(2)%record = "Why does the Universe bother to exist?"
        Decoration%List(3)%record = "What is the origin of mass and matter?"
        Decoration%List(4)%record = "What is the origin of life?"

        Decoration%text =   NLC//NLC// &
                            Decoration%List(1)%record // NLC//NLC// &
                            Decoration%List(2)%record // NLC//NLC// &
                            Decoration%List(3)%record // NLC//NLC// &
                            Decoration%List(4)%record // NLC//NLC//NLC

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"Test_Decoration_mod@test_writeDecoratedText_2."//num2str(Test%Image%id)//".out", status = "replace")

        call Decoration%writeDecoratedText  ( Decoration%text &
                                            , newLine = NLC &
                                            , outputUnit = fileUnit &
                                            )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)
            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            if (iostat/=0_IK) then
                assertion = .false.
                return
            end if
            OutputList(i)%record = trim(adjustl(OutputList(i)%record))

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if

        end do

    end function test_writeDecoratedText_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_writeDecoratedList_1() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical                         :: assertion
        logical                         :: assertionCurrent
        type(Decoration_type)           :: Decoration
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        integer(IK)                     :: fileUnit, i, iostat
        integer(IK), parameter          :: NLINE = 6_IK

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = "************************************************************************************************************************************"
        OutputList_ref(2)%record = "****                                                  Have you asked yourself:                                                  ****"
        OutputList_ref(3)%record = "****                                           Why does the Universe bother to exist?                                           ****"
        OutputList_ref(4)%record = "****                                           What is the origin of mass and matter?                                           ****"
        OutputList_ref(5)%record = "****                                                What is the origin of life?                                                 ****"
        OutputList_ref(6)%record = "************************************************************************************************************************************"

        allocate(Decoration%List(4))

        Decoration%List(1)%record = "Have you asked yourself:"
        Decoration%List(2)%record = "Why does the Universe bother to exist?"
        Decoration%List(3)%record = "What is the origin of mass and matter?"
        Decoration%List(4)%record = "What is the origin of life?"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"Test_Decoration_mod@test_writeDecoratedList_1."//num2str(Test%Image%id)//".out", status = "replace")

        call Decoration%writeDecoratedList(Decoration%List, outputUnit = fileUnit)

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)
            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            if (iostat/=0_IK) then
                assertion = .false.
                return
            end if
            OutputList(i)%record = trim(adjustl(OutputList(i)%record))

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if

        end do

    end function test_writeDecoratedList_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_writeDecoratedList_2() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical                         :: assertion
        logical                         :: assertionCurrent
        type(Decoration_type)           :: Decoration
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        integer(IK)                     :: fileUnit, i, iostat
        integer(IK), parameter          :: NLINE = 11_IK

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref( 1)%record = ""
        OutputList_ref( 2)%record = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        OutputList_ref( 3)%record = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        OutputList_ref( 4)%record = "%%                                                  Have you asked yourself:                                                  %%"
        OutputList_ref( 5)%record = "%%                                           Why does the Universe bother to exist?                                           %%"
        OutputList_ref( 6)%record = "%%                                           What is the origin of mass and matter?                                           %%"
        OutputList_ref( 7)%record = "%%                                                What is the origin of life?                                                 %%"
        OutputList_ref( 8)%record = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        OutputList_ref( 9)%record = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        OutputList_ref(10)%record = ""
        OutputList_ref(11)%record = ""

        allocate(Decoration%List(4))

        Decoration%List(1)%record = "Have you asked yourself:"
        Decoration%List(2)%record = "Why does the Universe bother to exist?"
        Decoration%List(3)%record = "What is the origin of mass and matter?"
        Decoration%List(4)%record = "What is the origin of life?"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"Test_Decoration_mod@test_writeDecoratedList_2."//num2str(Test%Image%id)//".out", status = "replace")

        call Decoration%writeDecoratedList  ( Decoration%List &
                                            , symbol = "%" &
                                            , width = 128_IK &
                                            , thicknessHorz = 2_IK &
                                            , thicknessVert = 2_IK &
                                            , marginTop = 1_IK &
                                            , marginBot = 2_IK &
                                            , outputUnit = fileUnit &
                                            )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)
            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            if (iostat/=0_IK) then
                assertion = .false.
                return
            end if
            OutputList(i)%record = trim(adjustl(OutputList(i)%record))

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if

        end do

    end function test_writeDecoratedList_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_drawLine_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: line_ref =   "HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!&
                                                    &HelloWorld!HelloWorld!HelloWorld!HelloWorld!Hello"
        character(:), allocatable   :: line
        line = drawLine(symbol = "HelloWorld!", width = 115_IK)
        assertion = line == line_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line_ref =", line_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line     =", line
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_drawLine_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_drawLine_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: line_ref =   "HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!&
                                                    &HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!HelloWorld!"
        character(:), allocatable   :: line
        line = drawLine(symbol = "HelloWorld!")
        assertion = line == line_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line_ref =", line_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line     =", line
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_drawLine_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_drawLine_3() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: line_ref =   "******************************************************************&
                                                    &******************************************************************"
        character(:), allocatable   :: line
        line = drawLine()
        assertion = line == line_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line_ref =", line_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "line     =", line
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_drawLine_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sandwich_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: sandwichedText_ref = "%                       The absence of evidence is not evidence for absence.                       %"
        character(:), allocatable   :: sandwichedText
        sandwichedText = sandwich   ( text = "The absence of evidence is not evidence for absence." &
                                    , symbol = "%" &
                                    , width = 100_IK &
                                    , thicknessHorz = 1_IK &
                                    )
        assertion = sandwichedText == sandwichedText_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText_ref   =", sandwichedText_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText       =", sandwichedText
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_sandwich_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sandwich_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: sandwichedText_ref = "%%%%                    The absence of evidence is not evidence for absence.                    %%%%"
        character(:), allocatable   :: sandwichedText
        sandwichedText = sandwich   ( text = "The absence of evidence is not evidence for absence." &
                                    , symbol = "%" &
                                    , width = 100_IK &
                                    )
        assertion = sandwichedText == sandwichedText_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText_ref   =", sandwichedText_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText       =", sandwichedText
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_sandwich_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sandwich_3() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: sandwichedText_ref = "%%%%                                    The absence of evidence is not evidence for absence.                                    %%%%"
        character(:), allocatable   :: sandwichedText
        sandwichedText = sandwich   ( text = "The absence of evidence is not evidence for absence." &
                                    , symbol = "%" &
                                    )
        assertion = sandwichedText == sandwichedText_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText_ref   =", sandwichedText_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText       =", sandwichedText
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_sandwich_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sandwich_4() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: sandwichedText_ref = "****                                    The absence of evidence is not evidence for absence.                                    ****"
        character(:), allocatable   :: sandwichedText
        sandwichedText = sandwich( text = "The absence of evidence is not evidence for absence." )
        assertion = sandwichedText == sandwichedText_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText_ref   =", sandwichedText_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText       =", sandwichedText
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_sandwich_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sandwich_5() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        character(*), parameter     :: sandwichedText_ref = "****                                                                                                                            ****"
        character(:), allocatable   :: sandwichedText
        sandwichedText = sandwich()
        assertion = sandwichedText == sandwichedText_ref
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15,:,' '))")
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText_ref   =", sandwichedText_ref
            write(Test%outputUnit,"(*(g0.15,:,' '))") "sandwichedText       =", sandwichedText
            write(Test%outputUnit,"(*(g0.15,:,' '))")
        end if
    end function test_sandwich_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_wrapText() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical                         :: assertion, assertionCurrent
        type(CharVec_type), allocatable :: ListOfLines_ref(:)
        type(CharVec_type), allocatable :: ListOfLines(:)
        integer(IK) , parameter         :: nline_ref = 6_IK
        character(:), allocatable       :: string
        integer(IK)                     :: nline, i

        assertion = .true.

        string =    "ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective &
                    &functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in &
                    &data science, Machine Learning, and scientific inference, with the design goal of unifying the &
                    &automation (of Monte Carlo simulations), user-friendliness (of the library), accessibility &
                    &(from multiple programming environments), high-performance (at runtime), and scalability &
                    &(across many parallel processors)."

        if (allocated(ListOfLines_ref)) deallocate(ListOfLines_ref); allocate(ListOfLines_ref(nline_ref))
        ListOfLines_ref(1)%record = "ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective "
        ListOfLines_ref(2)%record = "functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in "
        ListOfLines_ref(3)%record = "data science, Machine Learning, and scientific inference, with the design goal of unifying the "
        ListOfLines_ref(4)%record = "automation (of Monte Carlo simulations), user-friendliness (of the library), accessibility (from "
        ListOfLines_ref(5)%record = "multiple programming environments), high-performance (at runtime), and scalability (across many "
        ListOfLines_ref(6)%record = "parallel processors)."

        ListOfLines = wrapText(string = string, width = 100_IK, split = " ", pad = "    ")
        nline = size(ListOfLines)

        assertion = nline == nline_ref

        if (.not. assertion) then

            if (Test%isDebugMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "nline_ref = ", nline_ref
                write(Test%outputUnit,"(*(g0))") "nline     = ", nline
                write(Test%outputUnit,"(*(g0))")
            end if

            return

        else

            do i = 1, nline

                assertionCurrent = ListOfLines(i)%record == ListOfLines_ref(i)%record
                assertion = assertion .and. assertionCurrent

                if (Test%isDebugMode .and. .not. assertionCurrent) then
                    write(Test%outputUnit,"(*(g0))")
                    write(Test%outputUnit,"(*(g0))") "ListOfLines_ref(",num2str(i),")%record = '", ListOfLines_ref(i)%record, "'"
                    write(Test%outputUnit,"(*(g0))") "ListOfLines    (",num2str(i),")%record = '", ListOfLines(i)%record, "'"
                    write(Test%outputUnit,"(*(g0))")
                end if

            end do

        end if

    end function test_wrapText

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Decoration_mod
