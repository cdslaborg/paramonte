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

module Test_Batse_mod

    !use Constants_mod, only: IK, RK
    use Batse_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Batse

    type(Test_type) :: Test

    integer(IK) , parameter :: np = 86
    real(RK)    , parameter :: LOG10PBOL = 0._RK
    real(RK)    , parameter :: LOG10EPK_LOG10PH(2,np)   = reshape( &
                                                        [ -3.0000000000000000_RK, 4.92000000000000_RK &
                                                        , -2.9000000000000000_RK, 4.92271550938680_RK &
                                                        , -2.8000000000000000_RK, 4.94170723276800_RK &
                                                        , -2.7000000000000000_RK, 4.96223780646120_RK &
                                                        , -2.6000000000000000_RK, 4.98414439011840_RK &
                                                        , -2.5000000000000000_RK, 5.00727609375000_RK &
                                                        , -2.4000000000000000_RK, 5.03149354874880_RK &
                                                        , -2.3000000000000000_RK, 5.05666847891400_RK &
                                                        , -2.2000000000000000_RK, 5.08268327147520_RK &
                                                        , -2.1000000000000000_RK, 5.10943054811640_RK &
                                                        , -2.0000000000000000_RK, 5.13681273600000_RK &
                                                        , -1.9000000000000000_RK, 5.16474163879080_RK &
                                                        , -1.8000000000000000_RK, 5.19313800768000_RK &
                                                        , -1.7000000000000000_RK, 5.22193111240920_RK &
                                                        , -1.6000000000000000_RK, 5.25105831229440_RK &
                                                        , -1.5000000000000000_RK, 5.28046462725000_RK &
                                                        , -1.4000000000000000_RK, 5.31010230881280_RK &
                                                        , -1.3000000000000000_RK, 5.33993041116600_RK &
                                                        , -1.2000000000000000_RK, 5.36991436216320_RK &
                                                        , -1.1000000000000000_RK, 5.40002553435240_RK &
                                                        , -1.0000000000000000_RK, 5.43024081600000_RK &
                                                        , -0.9000000000000000_RK, 5.46054218211480_RK &
                                                        , -0.8000000000000000_RK, 5.49091626547200_RK &
                                                        , -0.7000000000000000_RK, 5.52135392763720_RK &
                                                        , -0.6000000000000000_RK, 5.55184982999040_RK &
                                                        , -0.5000000000000000_RK, 5.58240200475000_RK &
                                                        , -0.4000000000000000_RK, 5.61301142599680_RK &
                                                        , -0.3000000000000000_RK, 5.64368158069800_RK &
                                                        , -0.2000000000000000_RK, 5.67441803973120_RK &
                                                        , -0.0999999999999996_RK, 5.70522802890840_RK &
                                                        , 0.00000000000000000_RK, 5.73612000000000_RK &
                                                        , 0.10000000000000000_RK, 5.76710320175880_RK &
                                                        , 0.20000000000000000_RK, 5.79818725094400_RK &
                                                        , 0.30000000000000000_RK, 5.82938170334520_RK &
                                                        , 0.40000000000000000_RK, 5.86069562480640_RK &
                                                        , 0.50000000000000000_RK, 5.89213716225000_RK &
                                                        , 0.60000000000000000_RK, 5.92371311470080_RK &
                                                        , 0.70000000000000000_RK, 5.95542850431000_RK &
                                                        , 0.80000000000000000_RK, 5.98728614737920_RK &
                                                        , 0.90000000000000000_RK, 6.01928622538440_RK &
                                                        , 1.00000000000000000_RK, 6.05142585600000_RK &
                                                        , 1.10000000000000000_RK, 6.08369866412280_RK &
                                                        , 1.20000000000000000_RK, 6.11609435289600_RK &
                                                        , 1.30000000000000000_RK, 6.14859827473320_RK &
                                                        , 1.40000000000000000_RK, 6.18119100234240_RK &
                                                        , 1.50000000000000000_RK, 6.21425921875001_RK &
                                                        , 1.60000000000000000_RK, 6.24740828569600_RK &
                                                        , 1.70000000000000000_RK, 6.27826369177797_RK &
                                                        , 1.80000000000000000_RK, 6.30171825203203_RK &
                                                        , 1.90000000000000000_RK, 6.31516935429406_RK &
                                                        , 2.00000000000000000_RK, 6.31782000000002_RK &
                                                        , 2.10000000000000000_RK, 6.30989986450613_RK &
                                                        , 2.20000000000000000_RK, 6.29197575116811_RK &
                                                        , 2.30000000000000000_RK, 6.26452081342198_RK &
                                                        , 2.40000000000000000_RK, 6.22791191910419_RK &
                                                        , 2.50000000000000000_RK, 6.18562562500000_RK &
                                                        , 2.60000000000000000_RK, 6.13382209600000_RK &
                                                        , 2.70000000000000000_RK, 6.07721678100000_RK &
                                                        , 2.80000000000000000_RK, 6.01662553600000_RK &
                                                        , 2.90000000000000000_RK, 5.95282848100000_RK &
                                                        , 3.00000000000000000_RK, 5.88657000000000_RK &
                                                        , 3.10000000000000000_RK, 5.81855874100000_RK &
                                                        , 3.20000000000000000_RK, 5.74946761600000_RK &
                                                        , 3.30000000000000000_RK, 5.67993380100000_RK &
                                                        , 3.40000000000000000_RK, 5.61055873600000_RK &
                                                        , 3.50000000000000000_RK, 5.54190812500000_RK &
                                                        , 3.60000000000000000_RK, 5.47451193600000_RK &
                                                        , 3.70000000000000000_RK, 5.40886440100000_RK &
                                                        , 3.80000000000000000_RK, 5.34542401600000_RK &
                                                        , 3.90000000000000000_RK, 5.28461354100000_RK &
                                                        , 4.00000000000000000_RK, 5.22495000000000_RK &
                                                        , 4.10000000000000000_RK, 5.17119063610000_RK &
                                                        , 4.20000000000000000_RK, 5.12504318720001_RK &
                                                        , 4.30000000000000000_RK, 5.08588032630002_RK &
                                                        , 4.40000000000000000_RK, 5.05301458240003_RK &
                                                        , 4.50000000000000000_RK, 5.02572031250004_RK &
                                                        , 4.60000000000000000_RK, 5.00325567360003_RK &
                                                        , 4.70000000000000000_RK, 4.98488459470005_RK &
                                                        , 4.80000000000000000_RK, 4.96989874879999_RK &
                                                        , 4.90000000000000000_RK, 4.95763952489999_RK &
                                                        , 5.00000000000000000_RK, 4.94752000000006_RK &
                                                        , 5.10000000000000000_RK, 4.93904691110002_RK &
                                                        , 5.20000000000000000_RK, 4.93184262720004_RK &
                                                        , 5.30000000000000000_RK, 4.92566712130003_RK &
                                                        , 5.40000000000000000_RK, 4.92043994240004_RK &
                                                        , 5.50000000000000000_RK, 4.92000000000000_RK ],shape=shape(LOG10EPK_LOG10PH))

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************


    subroutine test_Batse()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_getLog10PF53()
        call test_getLogPF53()
        call Test%finalize()
        
    end subroutine test_Batse

!***********************************************************************************************************************************
!***********************************************************************************************************************************


    subroutine test_getLog10PF53()

        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)         :: ip
        character(len=1000) :: msg
        real(RK), parameter :: tolerance = 1.e-12_RK

        call Test%testing("getLog10PF53()")

        do ip = 1,np
            Test%assertion = abs(LOG10EPK_LOG10PH(2,ip)-getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL))<tolerance
            write(msg,* ) "The error with respect to reference value is larger than the tolerance.", new_line("a") &
                        , "tolerance, ip, LOG10EPK_LOG10PH(2,ip), getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL): " &
                        , tolerance, ip, LOG10EPK_LOG10PH(2,ip), getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL)
            call Test%verify(trim(adjustl(msg)))
        end do

    end subroutine test_getLog10PF53

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getLogPF53()

        use Constants_mod, only: IK, RK, LN10
        implicit none
        integer(IK) , parameter :: np = 86
        integer(IK)         :: ip
        character(len=1000) :: msg
        real(RK), parameter :: tolerance = 1.e-12_RK

        call Test%testing("getLogPF53()")

        do ip = 1,np
            Test%assertion = abs( getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL) &
                                -   getLogPF53(LOG10EPK_LOG10PH(1,ip)*LN10,LOG10PBOL)/LN10 &
                                ) < tolerance
            write(msg,* ) "The error with respect to reference value is larger than the tolerance.", new_line("a") &
                        , "tolerance, ip, LOG10EPK_LOG10PH(2,ip), getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL): " &
                        , tolerance, ip, LOG10EPK_LOG10PH(2,ip), getLog10PF53(LOG10EPK_LOG10PH(1,ip),LOG10PBOL)
            call Test%verify(trim(adjustl(msg)))
        end do

    end subroutine test_getLogPF53

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Batse_mod