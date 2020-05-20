!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module Test_Math_mod

    use Math_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Math

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Math()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_getLowerGamma()
        call test_getUpperGamma()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif
    end subroutine test_Math

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getLowerGamma()

        use Constants_mod, only: RK, IK

        implicit none
        integer, parameter      :: ntest = 3
        real(RK)                :: Exponent(ntest), UpperLim(ntest), LowerGamma(ntest)
        real(RK), parameter     :: LowerGammaRef(ntest) = [0.632120558828558_RK,0.184736755476228_RK,0.999817189367018_RK]
        real(RK)                :: difference
        integer                 :: i

        if (Test%Image%isFirst) call Test%testing("getLowerGamma")

        Exponent = [1.0_RK, 5.0_RK, 0.5_RK]
        UpperLim = [1.0_RK, 3.0_RK, 7.0_RK]

        !***************************************************************************************************************************

        do i = 1,ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            , tolerance = 1.e-7_RK &
                                            )
            difference = 2 * abs(LowerGammaRef(i) - LowerGamma(i)) / (LowerGammaRef(i) + LowerGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < 1.e-7_RK
            call Test%verify()
        end do
        !call Test%skipping()

        !***************************************************************************************************************************

        do i = 1,ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                            , tolerance = 1.e-3_RK &
                                            )
            difference = 2 * abs(LowerGammaRef(i) - LowerGamma(i)) / (LowerGammaRef(i) + LowerGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < 1.e-3_RK
            call Test%verify()
        end do
        !call Test%skipping()

        !***************************************************************************************************************************

        do i = 1,ntest
            LowerGamma(i) = getLowerGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , upperLim = UpperLim(i) &
                                           !, tolerance = 1.e-3_RK &
                                            )
            difference = 2 * abs(LowerGammaRef(i) - LowerGamma(i)) / (LowerGammaRef(i) + LowerGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, UpperLim, Reference LowerGamma, Computed LowerGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), UpperLim(i), LowerGamma(i), LowerGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < epsilon(1._RK) * 1.e1_RK
            call Test%verify()
        end do
        !call Test%skipping()

    end subroutine test_getLowerGamma

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getUpperGamma()

        use Constants_mod, only: RK, IK

        implicit none
        integer, parameter      :: ntest = 3
        real(RK)                :: Exponent(ntest), LowerLim(ntest), UpperGamma(ntest)
        real(RK), parameter     :: UpperGammaRef(ntest) = [0.367879441171442_RK,0.815263244523772_RK,1.828106329818355e-04_RK]
        real(RK)                :: difference
        integer                 :: i

        if (Test%Image%isFirst) call Test%testing("getUpperGamma")

        Exponent = [1.0_RK, 5.0_RK, 0.5_RK]
        LowerLim = [1.0_RK, 3.0_RK, 7.0_RK]

        !***************************************************************************************************************************

        do i = 1,ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , LowerLim = LowerLim(i) &
                                            , tolerance = 1.e-7_RK &
                                            )
            difference = 2 * abs(UpperGammaRef(i) - UpperGamma(i)) / (UpperGammaRef(i) + UpperGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, LowerLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < 1.e-7_RK
            call Test%verify()
        end do
        !call Test%skipping()

        !***************************************************************************************************************************

        do i = 1,ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , LowerLim = LowerLim(i) &
                                            , tolerance = 1.e-3_RK &
                                            )
            difference = 2 * abs(UpperGammaRef(i) - UpperGamma(i)) / (UpperGammaRef(i) + UpperGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, LowerLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < 1.e-3_RK
            call Test%verify()
        end do
        !call Test%skipping()

        !***************************************************************************************************************************

        do i = 1,ntest
            UpperGamma(i) = getUpperGamma   ( exponent = Exponent(i) &
                                            , logGammaExponent = log_gamma(Exponent(i)) &
                                            , LowerLim = LowerLim(i) &
                                           !, tolerance = 1.e-3_RK &
                                            )
            difference = 2 * abs(UpperGammaRef(i) - UpperGamma(i)) / (UpperGammaRef(i) + UpperGamma(i))
            if (Test%isDebugMode .and. Test%Image%isFirst) then
                write(Test%outputUnit,"(*(g0,:,', '))")
                write(Test%outputUnit,"(*(g0,:,', '))") "Exponent, LowerLim, Reference UpperGamma, Computed UpperGamma, difference:"
                write(Test%outputUnit,"(*(g0,:,', '))") Exponent(i), LowerLim(i), UpperGamma(i), UpperGammaRef(i), difference
                write(Test%outputUnit,"(*(g0,:,', '))")
            end if
            Test%assertion = difference < epsilon(1._RK) * 1.e2_RK
            call Test%verify()
        end do
        !call Test%skipping()

    end subroutine test_getUpperGamma

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Math_mod