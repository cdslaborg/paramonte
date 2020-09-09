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

module Test_CorrCoef_mod

    use CorrCoef_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_CorrCoef

    type(Test_type) :: Test
    integer, parameter  :: ndata = 50
    real(RK), parameter :: Data1(ndata) = &
                                        [ 120.131280000000_RK &
                                        , 118.789800000000_RK &
                                        , 118.895800000000_RK &
                                        , 119.379160000000_RK &
                                        , 117.751090000000_RK &
                                        , 115.396630000000_RK &
                                        , 113.495780000000_RK &
                                        , 114.635740000000_RK &
                                        , 117.231620000000_RK &
                                        , 118.123780000000_RK &
                                        , 119.710750000000_RK &
                                        , 117.265190000000_RK &
                                        , 114.134970000000_RK &
                                        , 121.689170000000_RK &
                                        , 115.801760000000_RK &
                                        , 117.407520000000_RK &
                                        , 118.966600000000_RK &
                                        , 121.113780000000_RK &
                                        , 115.724250000000_RK &
                                        , 118.406980000000_RK &
                                        , 117.778200000000_RK &
                                        , 115.731490000000_RK &
                                        , 119.184340000000_RK &
                                        , 117.686110000000_RK &
                                        , 119.619380000000_RK &
                                        , 115.007830000000_RK &
                                        , 118.544290000000_RK &
                                        , 119.010710000000_RK &
                                        , 117.242960000000_RK &
                                        , 116.145870000000_RK &
                                        , 117.908530000000_RK &
                                        , 118.999680000000_RK &
                                        , 113.657890000000_RK &
                                        , 114.844810000000_RK &
                                        , 117.666480000000_RK &
                                        , 114.232200000000_RK &
                                        , 114.647730000000_RK &
                                        , 120.284660000000_RK &
                                        , 115.660780000000_RK &
                                        , 117.098840000000_RK &
                                        , 114.627280000000_RK &
                                        , 117.065790000000_RK &
                                        , 116.938190000000_RK &
                                        , 121.164640000000_RK &
                                        , 116.831740000000_RK &
                                        , 118.905820000000_RK &
                                        , 117.057260000000_RK &
                                        , 116.203920000000_RK &
                                        , 116.729220000000_RK &
                                        , 118.009200000000_RK ]
    real(RK), parameter :: Data2(ndata) = &
                                        [ 5.28935260000000_RK &
                                        , 5.50145870000000_RK &
                                        , 5.89022390000000_RK &
                                        , 5.06549460000000_RK &
                                        , 5.62128260000000_RK &
                                        , 4.49246930000000_RK &
                                        , 3.54559920000000_RK &
                                        , 4.17171310000000_RK &
                                        , 5.34432780000000_RK &
                                        , 4.30855910000000_RK &
                                        , 6.12466330000000_RK &
                                        , 4.45103540000000_RK &
                                        , 4.08259680000000_RK &
                                        , 7.64761290000000_RK &
                                        , 6.53095480000000_RK &
                                        , 6.07550490000000_RK &
                                        , 7.32100850000000_RK &
                                        , 5.82501650000000_RK &
                                        , 4.19347540000000_RK &
                                        , 4.89687790000000_RK &
                                        , 5.61290890000000_RK &
                                        , 5.70994940000000_RK &
                                        , 5.00047920000000_RK &
                                        , 5.47741520000000_RK &
                                        , 4.99151560000000_RK &
                                        , 5.08172850000000_RK &
                                        , 5.98773500000000_RK &
                                        , 6.97849360000000_RK &
                                        , 6.91612860000000_RK &
                                        , 4.90595890000000_RK &
                                        , 5.71852950000000_RK &
                                        , 4.12146660000000_RK &
                                        , 5.51241440000000_RK &
                                        , 5.26293780000000_RK &
                                        , 5.14932990000000_RK &
                                        , 4.14738170000000_RK &
                                        , 5.55786790000000_RK &
                                        , 7.08800450000000_RK &
                                        , 6.08987380000000_RK &
                                        , 4.73697940000000_RK &
                                        , 3.80934450000000_RK &
                                        , 6.03942270000000_RK &
                                        , 5.96600840000000_RK &
                                        , 6.06674510000000_RK &
                                        , 5.84361600000000_RK &
                                        , 6.19013970000000_RK &
                                        , 4.43891700000000_RK &
                                        , 4.45833300000000_RK &
                                        , 5.47659170000000_RK &
                                        , 4.65761920000000_RK ]

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_CorrCoef()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_getCorrCoefSpearman()
        call Test%finalize()

    end subroutine test_CorrCoef

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getCorrCoefSpearman()
        
        use System_mod, only: sleep
        use Constants_mod, only: RK, IK
        use CorrCoef_mod

        implicit none
        type(CorrCoefSpearman_type) :: Spearman
        real(RK)                    :: refCorrCoef = 0.443361344537815_RK, refCorrPval = 0.00140031209338936_RK

        if (Test%Image%isFirst) call Test%testing("getCorrCoefSpearman")

        call Spearman%get   ( ndata             = ndata                     &
                            , Data1             = Data1                     &
                            , Data2             = Data2                     &
                            , rho               = Spearman%rho              &
                            , rhoProb           = Spearman%rhoProb          &
                            , dStarStar         = Spearman%dStarStar        &
                            , dStarStarSignif   = Spearman%dStarStarSignif  &
                            , dStarStarProb     = Spearman%dStarStarProb    &
                            , Err               = Spearman%Err              &
                            )
        if (Spearman%Err%occurred) then
            Test%assertion = .false.
            call Test%verify()
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Computed Spearman's statistic:"
            write(Test%outputUnit,"(*(g0))") "rho               = ", Spearman%rho
            write(Test%outputUnit,"(*(g0))") "rhoProb           = ", Spearman%rhoProb
            write(Test%outputUnit,"(*(g0))") "dStarStar         = ", Spearman%dStarStar
            write(Test%outputUnit,"(*(g0))") "dStarStarSignif   = ", Spearman%dStarStarSignif
            write(Test%outputUnit,"(*(g0))") "dStarStarProb     = ", Spearman%dStarStarProb
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Reference Spearman's statistic:"
            write(Test%outputUnit,"(*(g0))") "rho               = ", refCorrCoef
            write(Test%outputUnit,"(*(g0))") "rhoProb           = ", refCorrPval
            write(Test%outputUnit,"(*(g0))")
        end if

        Test%assertion = abs(Spearman%rho - refCorrCoef) / (Spearman%rho + refCorrCoef) < 1.e-10_RK
        call Test%verify()
        Test%assertion = abs(Spearman%rhoProb-refCorrPval) / (Spearman%rhoProb+refCorrPval) < 1.e-1_RK
        call Test%verify()
        !call Test%skipping()

    end subroutine test_getCorrCoefSpearman

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_CorrCoef_mod