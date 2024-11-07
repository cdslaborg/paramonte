!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This file contains procedure implementations of [pm_quadTest](@ref pm_quadTest).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_quadTest) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_val2str, only: getStr
    use pm_except, only: getInfNeg
    use pm_except, only: getInfPos
    use pm_quadPack, only: GK15, GK21
    use pm_quadPack, only: GK31, GK41
    use pm_quadPack, only: GK51, GK61
    use pm_quadPack, only: getQuadErr
    use pm_quadPack, only: isFailedQuad
    use pm_quadPack, only: weps, pinf, ninf
    use pm_quadPack, only: setNodeWeightGK
    use pm_arrayResize, only: setResized
    use pm_strASCII, only: getStrUpper
    use pm_io, only: display_type
    use pm_option, only: getOption

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_isFailedQuad_ENABLED 1
    ! \bug gfortran cannot recognize the procedure arguments without duplicating the full interface in the submodule.
    module subroutine test_isFailedQuad_RKH(disp, integrand, abstol, reltol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: test_isFailedQuad_RKH
#endif
        use pm_kind, only: RKG => RKH
        use pm_io, only: display_type
        type(display_type)      , intent(inout)         :: disp
        class(integrand_type)   , intent(in)            :: integrand
        real(RKG)               , intent(in), optional  :: abstol, reltol
#include "pm_quadTest@routines.inc.F90"
    end subroutine
#undef test_isFailedQuad_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getQuadErr_ENABLED 1
    ! \bug gfortran cannot recognize the procedure arguments without duplicating the full interface in the submodule.
    module subroutine test_getQuadErr_RKH(disp, integrand, atol, rtol, nintmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: test_getQuadErr_RKH
#endif
        use pm_kind, only: RKG => RKH
        use pm_io, only: display_type
        type(display_type)      , intent(inout)         :: disp
        class(integrand_type)   , intent(in)            :: integrand
        real(RKG)               , intent(in), optional  :: atol, rtol
        integer(IK)             , intent(in), optional  :: nintmax
#include "pm_quadTest@routines.inc.F90"
    end subroutine
#undef test_getQuadErr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int1_typer
        use pm_except, only: getInfNeg, getInfPos
        use pm_option, only: getOption
        use pm_kind, only: RKG => RKH
        self%lb = getOption(getInfNeg(0._RKG), lb)
        self%ub = getOption(getInfPos(0._RKG), ub)
        self%integral   = (atan(self%ub/2._RKG) * 2._RKG - atan(self%ub)) / 3._RKG & ! LCOV_EXCL_LINE
                        - (atan(self%lb/2._RKG) * 2._RKG - atan(self%lb)) / 3._RKG
        self%desc = "int1_type: an algebraic integrand of the form f(x) = x**2 / (x**2 + 1) / (x**2 + 4) for x in (lb, ub)"
    end procedure

    module procedure getInt1
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getInt1(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = x**2 / (x**2 + 1._RKG) / (x**2 + 4._RKG)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int2_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%a = getOption(1._RKG, a)
        self%b = getOption(1._RKG, b)
        self%ub = self%a / self%b
        self%lb = 0._RKG
        self%integral = 2._RKG * (sqrt(self%a - self%b * self%lb) - sqrt(self%a - self%b * self%ub)) / self%b
        self%desc = "int2_type: an algebraic integrand of the form f(x) = 1 / sqrt(a - b * x) for x in (0, a / b) with a > 0 and b > 0 with a singularity at the upper bound of integration"
        CHECK_ASSERTION(__LINE__, self%a > 0._RKG, SK_"@int2_typer(): The condition `self%a > 0._RKG` must hold. self%lb, self%a = "//getStr(self%a))
        CHECK_ASSERTION(__LINE__, self%b > 0._RKG, SK_"@int2_typer(): The condition `self%b > 0._RKG` must hold. self%lb, self%a = "//getStr(self%b))
    end procedure

    module procedure getInt2
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getInt2(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = 1._RKG / sqrt(self%a - self%b * x)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int3_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = 0._RKG
        self%ub = getOption(+1._RKG, ub)
        CHECK_ASSERTION(__LINE__, 0._RKG <= self%ub, SK_"@int3_typer(): The condition `0._RKG <= self%ub` must hold. self%ub = "//getStr(self%ub))
        self%integral = 2._RKG * sqrt(self%ub) * (log(self%ub) - 2._RKG)
        self%desc = "int3_type: an algebraic integrand of the form f(x) = log(x) / sqrt(x) for x in (0, ub) with a singularity at the lower limit of integration"
    end procedure

    module procedure getInt3
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getInt3(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = log(x) / sqrt(x)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int4_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = 0._RKG
        self%ub = 1._RKG
        self%integral = -0.189275187882093321180367135892330338053417661540147291526012234_RKG
        self%desc = "int4_type: an algebraic integrand of the form f(x) = log(x) / (1. + log(x)**2)**2 for x in (0, 1) with a singularity at the lower limit"
    end procedure

    module procedure getInt4
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getInt4(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = log(x) / (1._RKG + log(x)**2)**2
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int5_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        real(RKG), parameter    :: breakList(*) = [-sqrt(2._RKG), -1._RKG, 1._RKG, sqrt(2._RKG)]
        integer(IK)             :: i
        self%lb = getOption(0._RKG, lb)
        self%ub = getOption(3._RKG, ub)
        self%break = [real(RKG) ::]
        do i = 1, size(breakList)
            if (self%lb <= breakList(i) .and. breakList(i) <= self%ub) self%break = [self%break, breakList(i)]
        end do
        self%integral = getIntegral(self%ub) - getIntegral(self%lb)
        self%desc = "int5_type: an algebraic integrand of the form f(x) = x**3 log(abs((x**2 - 1) * (x**2 - 2))) for x in ("//getTTZ(getStr([self%lb, self%ub]))//SK_") with 4 possible singularities: [-sqrt(2.), -1., 1., sqrt(2.)]"
    contains
        pure function getIntegral(x) result(integral)
            real(RKG), intent(in)   :: x
            real(RKG)               :: integral
            integral = 4._RKG * real(log(cmplx(x**2 - 2._RKG, 0._RKG, RKG)), RKG) + real(log(cmplx(x**2 - 1._RKG, 0._RKG, RKG)), RKG)
            integral = 0.25_RKG * (x**4 * log(abs((x**2 - 1._RKG) * (x**2 - 2._RKG))) - integral - 3._RKG * x**2 - x**4)
        end function
    end procedure

    module procedure getInt5
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getInt5(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        CHECK_ASSERTION(__LINE__, all(10 * abs(x - nearest(x, x)) < abs(x - self%break)), \
        SK_"@getInt5(): The condition `all(10 * abs(x - nearest(x, x)) < abs(x - self%break))` must hold. 10 * abs(x - nearest(x, x)), abs(x - self%break) = "//getStr([real(RKG) :: 10 * abs(x - nearest(x, x)), abs(x - self%break)]))
        func = x**3 * log(abs((x**2 - 1._RKG) * (x**2 - 2._RKG)))
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int6_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = 0._RKG
        self%ub = getInfPos(0._RKG)
        self%integral = -acos(-1._RKG) * log(10._RKG) / 20._RKG
        self%desc = "int6_type: an algebraic integrand of the form f(x) = log(x) / (1 + 100 * x**2) for x in (0, +Inf) with a singularity at the lower limit"
    end procedure

    module procedure getInt6
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getInt6(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = log(x) / (1._RKG + 100._RKG * x**2)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int7_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = 1._RKG / 3._RKG
        self%ub = getInfPos(0._RKG)
        self%break = [1._RKG / sqrt(2._RKG), 1._RKG]
        self%integral = 52.7407483834714449977291997202299809_RKG
        self%desc = "int7_type: an algebraic integrand of the form f(x) = -log(abs((1 - x**2) * (1 - 2 * x**2)) / x**4) / x**5 for x in (1./3., +Inf) with two singularities at [1 / sqrt(2), 1]"
    end procedure

    module procedure getInt7
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getInt7(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        CHECK_ASSERTION(__LINE__, all(10 * abs(x - nearest(x, x)) < abs(x - self%break)), \
        SK_"@getInt7(): The condition `all(10 * abs(x - nearest(x, x)) < abs(x - self%break))` must hold. 10 * abs(x - nearest(x, x)), abs(x - self%break) = "//getStr([real(RKG) :: 10 * abs(x - nearest(x, x)), abs(x - self%break)]))
        func = log(abs((1._RKG - x**2) * (1._RKG - 2 * x**2)) / x**4) / x**5
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int8_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getInfNeg(0._RKG)
        self%ub = -1._RKG / 3._RKG
        self%break = [-1._RKG, -1._RKG / sqrt(2._RKG)]
        self%integral = 52.7407483834714449977291997202299809_RKG
        self%desc = "int8_type: an algebraic integrand of the form f(x) = log(abs((1 - x**2) * (1 - 2 * x**2)) / x**4) / x**5 for x in (-Inf, -1./3.) with two singularities at [-1, -1 / sqrt(2)]"
    end procedure

    module procedure getInt8
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getInt8(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        CHECK_ASSERTION(__LINE__, all(10 * abs(x - nearest(x, x)) < abs(x - self%break)), \
        SK_"@getInt8(): The condition `all(10 * abs(x - nearest(x, x)) < abs(x - self%break))` must hold. 10 * abs(x - nearest(x, x)), abs(x - self%break) = "//getStr([real(RKG) :: 10 * abs(x - nearest(x, x)), abs(x - self%break)]))
        func = -log(abs((1._RKG - x**2) * (1._RKG - 2 * x**2)) / x**4) / x**5
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure int9_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getInfNeg(0._RKG)
        self%ub = getInfPos(0._RKG)
        self%break = [-10._RKG, -9._RKG, 1._RKG / 3._RKG, 1._RKG / sqrt(2._RKG), 1._RKG]
        self%integral = 53.7407483834714449977291997202299809_RKG
        self%desc = "int9_type: an algebraic piecewise integrand of the form f(x) = log(abs((1 - x**2) * (1 - 2 * x**2)) / x**4) / x**5 for x in (1./3., +Inf) and 1 / (acos(-1) * sqrt(-(x+10) * (x+9))) for x in (-10, 9), otherwise 0, with two singularities at [-10, -9, 1/3., 1 / sqrt(2), 1]"
    end procedure

    module procedure getInt9
        use pm_kind, only: RKG => RKH
        !check_assertion(__LINE__, all(10 * abs(x - nearest(x, x)) < abs(x - self%break)), \
        !SK_"@getInt9(): The condition `all(abs(x - nearest(x, x)) < abs(x - self%break))` must hold. abs(x - nearest(x, x)), abs(x - self%break) = "//getStr([real(RKG) :: abs(x - nearest(x, x)), abs(x - self%break)]))
        if (all(abs(x - nearest(x, x)) < 10 * abs(x - self%break))) then
            func = 0._RKG
            if (x > 1._RKG / 3._RKG) func = log(abs((1._RKG - x**2) * (1._RKG - 2 * x**2)) / x**4) / x**5
            if (-10._RKG < x .and. x < -9._RKG) func = 1._RKG / (acos(-1._RKG) * sqrt(-(x + 10._RKG) * (x + 9._RKG)))
        else
            func = 1.e9_RKG
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intGamUpp_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos, setNAN, isInfPos
        use pm_mathGamma, only: getGammaIncUpp
        real(RKG) :: lbint, ubint
        if (present(lb)) then
            CHECK_ASSERTION(__LINE__, 0 < lb, SK_"@getIntGamUpp(): The condition `0 < lb` must hold. lb = "//getStr(lb))
            self%lb = lb
        else
            self%lb = 1._RKG
        end if
        if (present(ub)) then
            CHECK_ASSERTION(__LINE__, 0 < ub, SK_"@getIntGamUpp(): The condition `0 < lb` must hold. lb = "//getStr(lb))
            self%ub = ub
        else
            self%ub = getInfPos(0._RKG)
        end if
        CHECK_ASSERTION(__LINE__, self%lb < self%ub, SK_"@getIntGamUpp(): The condition `lb < ub` must hold. lb, ub = "//getStr([self%lb, self%ub]))
        if (present(alpha)) then
            self%alpha = alpha
        else
            self%alpha = 1._RKG
        end if
        if (present(beta)) then
            CHECK_ASSERTION(__LINE__, 0 < beta, SK_"@getIntGamUpp(): The condition `0 < beta` must hold. beta = "//getStr(beta))
            self%beta = beta
        else
            self%beta = 1._RKG
        end if
        if (self%alpha > -1) then
            lbint = gamma(self%alpha + 1) * getGammaIncUpp(self%beta * self%lb, kappa = self%alpha + 1)
            lbint = lbint * exp(self%beta * self%lb) / (self%lb**self%alpha * self%beta**(self%alpha + 1))
            if (isInfPos(self%ub)) then
                ubint = 0._RKG
            else
                ubint = gamma(self%alpha + 1) * getGammaIncUpp(self%beta * self%ub, kappa = self%alpha + 1)
                ubint = ubint * exp(self%beta * self%ub) / (self%ub**self%alpha * self%beta**(self%alpha + 1))
            end if
            self%integral = lbint + ubint
        else
            call setNAN(self%integral)
        end if
        self%desc = "intGamUpp_type: an algebraic integrand of the form f(x; lb, alpha, beta) = (x / lb)**alpha * exp(-beta * (x - lb)) for x in (lb, +Inf), lb > 0."
    end procedure

    module procedure getIntGamUpp
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getIntGamUpp(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = (x / self%lb)**self%alpha * exp(-self%beta * (x - self%lb))
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intSinCos_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getOption(-1_IK, lf) * acos(-1._RKG)
        self%ub = getOption(+1_IK, uf) * acos(-1._RKG)
        self%a = getOption(1._RKG, a)
        self%b = getOption(1._RKG, b)
        self%integral = (self%ub - self%lb) * bessel_j0(self%a)
        self%desc = "intSinCos_typer(): a highly oscillatory integrand of the form f(x) = cos(a * sin(b * x)) for x in (lb, ub)"
    end procedure

    module procedure getIntSinCos
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntSinCos(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = cos(self%a * sin(self%b  * x))
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intNormPDF_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_distNorm, only: getNormCDF
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getOption(getInfNeg(0._RKG), lb)
        self%ub = getOption(getInfPos(0._RKG), ub)
        self%mu = getOption(+0._RKG, mu)
        self%sigma = getOption(1._RKG, sigma)
        self%invSigma = 1._RKG / self%sigma
        self%logInvSigma = log(self%invSigma)
        self%integral = getNormCDF(self%ub, self%mu, self%sigma) - getNormCDF(self%lb, self%mu, self%sigma)
        self%desc = "intNormPDF_type: f(x) = exp(-0.5 * (log(x) - mu)**2 / sigma**2) / (sigma * sqrt(2 * acos(-1.))) for x in (lb, ub)"
    end procedure

    module procedure getIntNormPDF
        use pm_kind, only: RKG => RKH
        use pm_distNorm, only: setNormLogPDF
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntNormPDF(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        call setNormLogPDF(func, x, self%mu, self%invSigma, self%logInvSigma)
        func = exp(func)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intLogNormPDF_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_distLogNorm, only: getLogNormCDF
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getOption(0._RKG, lb)
        self%ub = getOption(getInfPos(0._RKG), ub)
        CHECK_ASSERTION(__LINE__, 0._RKG <= self%lb .and. self%lb <= self%ub, \
        SK_"@intLogNormPDF_typer(): The condition `0._RKG <= self%lb .and. slef%lb <= self%ub` must hold. self%lb, self%ub = "//getStr([self%lb, self%ub]))
        self%mu = getOption(+0._RKG, mu)
        self%sigma = getOption(1._RKG, sigma)
        self%invSigma = 1._RKG / self%sigma
        self%logInvSigma = log(self%invSigma)
        self%integral = getLogNormCDF(self%ub, self%mu, self%sigma) - getLogNormCDF(self%lb, self%mu, self%sigma)
        self%desc = "intLogNormPDF_type: f(x) = exp(-0.5 * (log(x) - mu)**2 / sigma**2) / (x * sigma * sqrt(2 * acos(-1.))) for x in (0 <= lb, ub)"
    end procedure

    module procedure getIntLogNormPDF
        use pm_kind, only: RKG => RKH
        use pm_distLogNorm, only: setLogNormLogPDF
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntLogNormPDF(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        call setLogNormLogPDF(func, log(x), self%mu, self%invSigma, self%logInvSigma)
        func = exp(func)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intGenExpGammaPDF_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
       !use pm_distGenExpGamma, only: getGenExpGammaCDF
        use pm_distGenExpGamma, only: getGenExpGammaLogPDFNF
        self%lb = getOption(getInfNeg(0._RKG), lb)
        self%ub = getOption(getInfPos(0._RKG), ub)
        self%kappa = getOption(1._RKG, kappa)
        self%invOmega = getOption(1._RKG, invOmega)
        self%logSigma = getOption(0._RKG, logSigma)
        self%logPDFNF = getGenExpGammaLogPDFNF(self%kappa, self%invOmega)
        CHECK_ASSERTION(__LINE__, 0._RKG < self%kappa, SK_"@intGenExpGammaPDF_typer(): The condition `0._RKG < kappa` must hold. kappa = "//getStr(self%kappa))
        CHECK_ASSERTION(__LINE__, 0._RKG < self%invOmega, SK_"@intGenExpGammaPDF_typer(): The condition `0._RKG < invOmega` must hold. invOmega = "//getStr(self%invOmega))
        self%integral = 1._RKG !getGenExpGammaCDF(exp(self%ub), self%mu, self%sigma) - getGenExpGammaCDF(self%lb, self%mu, self%sigma)
        self%desc = SK_"intGenExpGammaPDF_type: f(x) = GenExpGamma(x; kappa = "//getTTZ(getStr(self%kappa))//SK_", omega = "//getTTZ(getStr(1._RKG/self%invOmega))//SK_", logSigma = "//getTTZ(getStr(self%logSigma))//SK_") for x in ("//getTTZ(getStr(self%lb))//SK_", "//getTTZ(getStr(self%ub))//SK_")"
    end procedure

    module procedure getIntGenExpGammaPDF
        use pm_kind, only: RKG => RKH
        use pm_distGenExpGamma, only: setGenExpGammaLogPDF
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getIntGenExpGammaPDF(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        call setGenExpGammaLogPDF(func, x, self%logPDFNF, self%kappa, self%invOmega, self%logSigma)
        func = exp(func)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intPentaGammaInf_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getInfNeg(0._RKG)
        self%ub = getInfPos(0._RKG)
        self%break = [-9._RKG, -5._RKG, 2._RKG, 5._RKG, 7._RKG]
        self%integral = 5._RKG
        self%desc = "intPentaGammaInf_type: f(x) = sum of five Gamma PDFs with five break points in the integration range x in (-Inf, +Inf)"
    end procedure

    module procedure getIntPentaGammaInf
        use pm_distGamma, only: getGammaLogPDF
        use pm_kind, only: RKG => RKH
        !check_assertion(__LINE__, all(10 * abs(x - nearest(x, x)) < abs(x - self%break)), \
        !SK_"@getIntPentaGammaInf(): The condition `all(10 * abs(x - nearest(x, x)) < abs(x - self%break))` must hold. 10 * abs(x - nearest(x, x)), abs(x - self%break) = "//getStr([real(RKG) :: 10 * abs(x - nearest(x, x)), abs(x - self%break)]))
        if (all(abs(x - nearest(x, x)) < 10 * abs(x - self%break))) then
            func = 0._RKG
            if (x > -9._RKG) func = func + exp(getGammaLogPDF(x + 9._RKG, 0.7_RKG, 1._RKG))
            if (x > -5._RKG) func = func + exp(getGammaLogPDF(x + 5._RKG, 0.7_RKG, 1._RKG))
            if (x > +5._RKG) func = func + exp(getGammaLogPDF(x - 5._RKG, 0.7_RKG, 1._RKG))
            if (x < +2._RKG) func = func + exp(getGammaLogPDF(2._RKG - x, 0.7_RKG, 1._RKG))
            if (x < +7._RKG) func = func + exp(getGammaLogPDF(7._RKG - x, 0.7_RKG, 1._RKG))
        else
            func = 1.e9_RKG
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intDoncker1_typer
        use pm_kind, only: RKG => RKH
        use pm_except, only: getInfNeg, getInfPos
        use pm_option, only: getOption
        self%lb = getOption(0._RKG, lb)
        self%ub = getOption(getInfPos(0._RKG), ub)
        CHECK_ASSERTION(__LINE__, 0._RKG <= self%lb, SK_"@intDoncker1_typer(): The condition `0._RKG <= self%lb` must hold. self%ub = "//getStr(self%ub))
        CHECK_ASSERTION(__LINE__, self%lb < self%ub, SK_"@intDoncker1_typer(): The condition `self%lb <= self%ub` must hold. self%ub = "//getStr([self%lb, self%ub]))
        self%integral = 2._RKG * (atan(sqrt(self%ub)) - atan(sqrt(self%lb)))
        self%desc = "intDoncker1_type: f(x) = 1 / (1 + x) / sqrt(x) for x in (0 <= lb, ub) with a square-root singularity at 0"
    end procedure

    module procedure getIntDoncker1
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntDoncker1(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = 1._RKG / (1._RKG + x) / sqrt(x)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intDoncker2_typer
        use pm_kind, only: RKG => RKH
        use pm_option, only: getOption
        use pm_except, only: getInfNeg, getInfPos
        self%lb = getOption(getInfNeg(0._RKG), lb)
        self%ub = getOption(0._RKG, ub)
        CHECK_ASSERTION(__LINE__, self%ub <= 0._RKG, SK_"@intDoncker2_typer(): The condition `0._RKG <= self%ub` must hold. self%ub = "//getStr(self%ub))
        CHECK_ASSERTION(__LINE__, self%lb < self%ub, SK_"@intDoncker2_typer(): The condition `self%lb <= self%ub` must hold. self%ub = "//getStr([self%lb, self%ub]))
        self%integral = sqrt(acos(-1._RKG)) * (erf(sqrt(-self%lb)) - erf(sqrt(-self%ub)))
        self%desc = "intDoncker2_type: f(x) = exp(x) / sqrt(-x) for x in (lb, ub <= 0) with a square-root singularity at 0"
    end procedure

    module procedure getIntDoncker2
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getIntDoncker2(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = exp(x) / sqrt(-x)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intCauchy1_typer
        use pm_kind, only: RKG => RKH
        self%lb = getOption(-2._RKG, lb)
        self%ub = getOption(+3._RKG, ub)
        self%wcauchy = wcauchy_type(getOption(+1._RKG, cs))
        CHECK_ASSERTION(__LINE__, self%lb < self%wcauchy%cs .and. self%wcauchy%cs < self%ub, SK_"@intCauchy1_typer(): The condition `self%lb < self%wcauchy%cs .and. self%wcauchy%cs < self%ub` must hold. self%lb, self%wcauchy%cs, self%ub = "//getStr([self%lb, self%wcauchy%cs, self%ub]))
        self%integral = log(self%ub - self%wcauchy%cs) - log(self%wcauchy%cs - self%lb)
        self%desc = "intCauchy1_type: an integrand of the form w(x) * f(x) with Cauchy weight w(x) 1 / (x - cs) and f(x) = 1 ~,~ x \in (lb < cs, cs < ub)"
    end procedure

    module procedure getIntCauchy1
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntCauchy1(): The condition `self%lb < x .and. x < self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        CHECK_ASSERTION(__LINE__, abs(x - nearest(x, x)) < abs(x - self%wcauchy%cs), \
        SK_"@getIntCauchy1(): The condition `abs(x - nearest(x, x)) < abs(x - self%wcauchy%cs)` must hold. 10 * abs(x - nearest(x, x)), abs(x - self%wcauchy%cs) = "//getStr([abs(x - nearest(x, x)), abs(x - self%wcauchy%cs)]))
        func = 1._RKG
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure intCauchy2_typer
        use pm_kind, only: RKG => RKH
        use pm_val2str, only: getStr
        use pm_swap, only: setSwapped
        use pm_mathMinMax, only: setMinMax
        self%lb = getOption(-3._RKG, lb)
        self%ub = getOption(+2._RKG, ub)
        self%pole(1) = getOption(-2._RKG, cs1)
        self%pole(2) = getOption(+3._RKG, cs2)
        call setMinMax(self%pole)
        CHECK_ASSERTION(__LINE__, .not. (self%lb < self%pole(1) .and. self%pole(1) < self%ub .and. self%lb < self%pole(2) .and. self%pole(2) < self%ub), \
        SK_"@getIntCauchy2(): The two specified Cauchy poles of the integrand must not simultaneously appear with the integration range. lb, cs1, cs2, ub = "//\
        getStr([self%lb, self%pole, self%ub]))
        self%wcauchy = wcauchy_type(self%pole(1))
        self%csnot = self%pole(2)
        if (self%lb < self%csnot .and. self%csnot < self%ub) call setSwapped(self%wcauchy%cs, self%csnot)
        self%integral = getIntegralIndef(self%ub) - getIntegralIndef(self%lb)
        self%desc = SK_"intCauchy2_type: an integrand of the form w(x) * f(x) = 1 / (x "//merge(SK_"+ ", SK_"- ", self%pole(1) < 0._RKG)//getTTZ(getStr(abs(self%pole(1))))// & ! LCOV_EXCL_LINE
        SK_") / (x "//merge(SK_"+ ", SK_"- ", self%pole(2) < 0._RKG)//getTTZ(getStr(abs(self%pole(2))))// & ! LCOV_EXCL_LINE
        SK_") for x in ("//getTTZ(getStr(self%lb))//SK_", "//getTTZ(getStr(self%ub))//SK_")"
    contains
        PURE function getIntegralIndef(x) result(integralIndef)
            real(RKG), intent(in)   :: x
            real(RKG)               :: integralIndef
            !integralIndef = real((log(cmplx(x - self%pole(1), 0._RKG, RKG)) - log(cmplx(x - self%pole(2), 0._RKG, RKG))), RKG) / (self%pole(1) - self%pole(2))
            integralIndef = real((log(cmplx(1._RKG + (self%pole(2) - self%pole(1)) / (x - self%pole(2)), 0._RKG, RKG))), RKG) / (self%pole(1) - self%pole(2))
        end function
    end procedure

    module procedure getIntCauchy2
        use pm_kind, only: RKG => RKH
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, \
        SK_"@getIntCauchy2(): The condition `self%lb < x .and. x < self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        !check_assertion(__LINE__, all(abs(x - nearest(x, x)) < abs(x - self%pole)), \
        !SK_"@getIntCauchy2(): The condition `all(10 * abs(x - nearest(x, x)) < abs(x - self%pole))` must hold. abs(x - nearest(x, x)), abs(x - self%pole)) = "//getStr([real(RKG) :: abs(x - nearest(x, x)), abs(x - self%pole)]))
        if (all(abs(x - nearest(x, x)) < abs(x - self%pole))) then
            func = 1._RKG / (x - self%csnot)
        else
            func = 1.e9_RKG
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines