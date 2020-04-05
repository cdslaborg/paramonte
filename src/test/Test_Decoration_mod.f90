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
