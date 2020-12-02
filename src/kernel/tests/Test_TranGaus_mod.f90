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

!>  \brief This module contains tests of the module [TranGaus_mod](@ref trangaus_mod).
!>  @author Amir Shahmoradi

module Test_TranGaus_mod

    use TranGaus_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_TranGaus

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_TranGaus()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getTranGaus_1, "test_getTranGaus_1")
        call Test%finalize()

    end subroutine test_TranGaus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getTranGaus_1() result(assertion)

        use Constants_mod, only: RK, IK
        use String_mod, only: num2str

        implicit none
        logical                 :: assertion
        integer(IK), parameter  :: NSIM = 100_IK
        real(RK), parameter     :: avg = 2._RK
        real(RK), parameter     :: std = 3._RK
        integer(IK)             :: isim, fileUnit1, fileUnit2, fileUnit3, fileUnit4, fileUnit5

        assertion = .true.

        open(newunit=fileUnit1, file=Test%outDir//"Test_TranGaus_mod@test_getTranGaus_Low1Upp9Avg2Std3."//num2str(Test%Image%id)//".temp", status="replace")
        open(newunit=fileUnit2, file=Test%outDir//"Test_TranGaus_mod@test_getTranGaus_Low20Upp30Avg2Std3."//num2str(Test%Image%id)//".temp", status="replace")
        open(newunit=fileUnit3, file=Test%outDir//"Test_TranGaus_mod@test_getTranGaus_Low-15Upp4Avg2Std3."//num2str(Test%Image%id)//".temp", status="replace")
        open(newunit=fileUnit4, file=Test%outDir//"Test_TranGaus_mod@test_getTranGaus_Low50Upp55Avg2Std3."//num2str(Test%Image%id)//".temp", status="replace")
        open(newunit=fileUnit5, file=Test%outDir//"Test_TranGaus_mod@test_getTranGaus_Low-50Upp55Avg2Std3."//num2str(Test%Image%id)//".temp", status="replace") ! test `lowerLim < xmin`
        write(fileUnit1,*) "TranGaus"
        write(fileUnit2,*) "TranGaus"
        write(fileUnit3,*) "TranGaus"
        write(fileUnit4,*) "TranGaus"

        ! Random variable generation

        do isim = 1, NSIM
            write(fileUnit1,*) getTranGaus(lowerLim=1._RK,upperLim=9._RK,avg=avg,std=std)
            write(fileUnit2,*) getTranGaus(lowerLim=20._RK,upperLim=30._RK,avg=avg,std=std)
            write(fileUnit3,*) getTranGaus(lowerLim=-15._RK,upperLim=4._RK,avg=avg,std=std)
            write(fileUnit4,*) getTranGaus(lowerLim=50._RK,upperLim=55._RK,avg=avg,std=std)
            write(fileUnit5,*) getTranGaus(lowerLim=-50._RK,upperLim=55._RK,avg=avg,std=std)
        end do
        close(fileUnit1)
        close(fileUnit2)
        close(fileUnit3)
        close(fileUnit4)
        close(fileUnit5)

        ! ideally, as a test, one could compute the PDF of the Gaussian and compare it to the histogram of the generated points.
        ! This is not an easy task. For now, we rust the giants.

    end function test_getTranGaus_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_TranGaus_mod ! LCOV_EXCL_LINE