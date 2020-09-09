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

module Test_Decoration_mod

    !use, intrinsic :: iso_fortran_env, only: output_unit
    use Test_mod, only: Test_type
    use Constants_mod, only: IK
    use Decoration_mod

    implicit none

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Decoration()
        implicit none
        type(Decoration_type) :: Decor

        Test = Test_type(moduleName=MODULE_NAME)

        if (Test%Image%isFirst) call Test%testing("Decoration_type")

        if (Test%isDebugMode .and. Test%Image%isFirst) then

            allocate(Decor%List(4))
            Decor%List(1)%record = "Have you asked yourself:"
            Decor%List(2)%record = "Why does the Universe bother to exist?"
            Decor%List(3)%record = "What is the origin of mass and matter?"
            Decor%List(4)%record = "What is the origin of life?"
            call Decor%writeDecoratedList(Decor%List)

            Decor%text = "\n\n" // &
                         Decor%List(1)%record // "\n\n" // &
                         Decor%List(2)%record // "\n\n" // &
                         Decor%List(3)%record // "\n\n" // &
                         Decor%List(4)%record // "\n\n\n"
            call Decor%writeDecoratedText( Decor%text &
                                         , newLine="\n" &
                                         , width = 32 &
                                         , symbol="&" &
                                         , thicknessHorz=4 &
                                         , thicknessVert=2 &
                                         , marginTop=2 &
                                         , marginBot=1 &
                                         )

            call write(Test%outputUnit,1,0,1, "Here is the output of drawLine('HelloWorld!'):" )
            call write(Test%outputUnit,1,0,1, drawLine(symbol="HelloWorld!",width=132) )

            call write(Test%outputUnit,1,0,1, "Here is the output of sandwich():" )
            call write(Test%outputUnit,1,0,1, sandwich() )

            call write(Test%outputUnit,1,0,1, "Here is the output of sandwich( text='Hello World!', symbol='@', width=20 ) :" )
            call write(Test%outputUnit,1,0,1, sandwich(text='Hello World!',symbol='@',width=20) )

            call write(Test%outputUnit,1,0,1, "Here is the output of sandwich( text='Hello World!', symbol='@', width=10 ) :" )
            call write(Test%outputUnit,1,0,1, sandwich(text='Hello World!',symbol='@',width=10) )

            call write(Test%outputUnit,1,0,1, "Here is the output of sandwich( text='Hello World!', symbol='@', width=10 ) :" )
            call write(Test%outputUnit,1,0,1, sandwich(text='Hello World!',symbol='@',width=5) )

            write(Test%outputUnit,"(A)") 

        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()
        call Test%finalize()


    end subroutine test_Decoration

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Decoration_mod
