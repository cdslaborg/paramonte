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
!>  This include file contains implementations of the procedures in module [pm_mathGammaGil](@ref pm_mathGammaGil).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, July 22, 2024, 11:45 AM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! incomplete gamma dependencies:
! qfraction
! qtaylor
! ptaylor
! pqasymp
! dompart
! alfa
!
! inverse incomplete gamma dependencies:
! lambdaeta
! sqrttwopi
! gamstar
! loggam
! inverfc
! eps1
! eps2
! eps3
! incgam

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGammaIncLowGil_ENABLED || getGammaIncUppGil_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     getGammaIncLowGil_ENABLED
        real(RKG) :: gammaIncUpp
        character(*, SK), parameter :: PROCEDURE_NAME = "@getGammaIncLowGil()"
#elif   getGammaIncUppGil_ENABLED
        real(RKG) :: gammaIncLow
        character(*, SK), parameter :: PROCEDURE_NAME = "@getGammaIncUppGil()"
#else
#error  "Unrecognized interface."
#endif
        call setGammaIncGil(gammaIncLow, gammaIncUpp, x, kappa, info)
        if (info < 0_IK) error stop SK_"@file::"//__FILE__//SK_"@line::"//getStr(__LINE__)//&! LCOV_EXCL_LINE
        MODULE_NAME//PROCEDURE_NAME//SK_": The Incomplete Gamma function failed to converge." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaIncGil_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: PI = acos(-1._RKG)
        real(RKG), parameter :: SQRT_2PI = sqrt(2._RKG * PI)
        real(RKG), parameter :: LOG_SQRT_2PI = log(SQRT_2PI)
        real(RKG), parameter :: EPS = 10 * epsilon(x) !< The default relative accuracy.
        real(RKG), parameter :: FPMAX = huge(x) * sqrt(EPS) !< A number near the largest representable floating-point number.
        real(RKG), parameter :: FPMIN = tiny(x) / sqrt(EPS) !< A number near the smallest representable floating-point number.
        real(RKG), parameter :: LOG_FPMIN = log(FPMIN) !< A number near the smallest representable floating-point number.
        real(RKG) :: logx, dp

        if (x < FPMIN) then
          logx = LOG_FPMIN
        else
          logx = log(x)
        end if

        info = 0_IK
        if (alfa(x) < kappa) then

            dp = dompart(kappa, x, .false._IK)
            if (dp < 0._RKG) then
                info = -1_IK
            else
                if ((x < 0.3_RKG * kappa) .or. (kappa < 12._RKG)) then
                    gammaIncLow = ptaylor(kappa, x, dp)
                else
                    gammaIncLow = pqasymp(kappa, x, dp, .true._LK)
                end if
                gammaIncUpp = 1._RKG - gammaIncLow
            end if

        else

            if (kappa < -FPMIN / logx) then
                gammaIncLow = 1._RKG
                gammaIncUpp = 0._RKG
            else
                if (x < 1._RKG) then
                    dp = dompart(kappa, x, .true._LK)
                    if (dp < 0._RKG) then
                        info = -1_IK
                    else
                        gammaIncUpp = qtaylor(kappa, x, dp)
                        gammaIncLow = 1._RKG - gammaIncUpp
                    end if
                else
                    dp = dompart(kappa, x, .false._LK)
                    if (dp < 0._RKG) then
                        info = -1_IK
                    else
                        if (x > 2.35_RKG * kappa .or. kappa < 12._RKG) then
                            gammaIncUpp = qfraction(kappa, x, dp)
                        else
                            gammaIncUpp = pqasymp(kappa, x, dp, .false._LK)
                        end if
                        gammaIncLow = 1._RKG - gammaIncUpp
                    end if
                end if

            end if

        end if

    contains

        pure function lnec(x) result(y)
            real(RKG), intent(in) :: x
            real(RKG) :: y, z, e2, r, s
            ! x > -1: y := ln1 := ln(1 + x) - x
            z = logoneplusx(x)
            y = z - x
            e2 = exmin1minx(z)
            s = .5_RKG * e2 * z * z
            r = (s + y) / (s + 1._RKG + z)
            y = y - r * (6._RKG - r) / (6._RKG - 4._RKG * r)
            ! y := x + y
        end function lnec

        pure function alfa(x) result(y)
            real(RKG), intent(in) :: x
            real(RKG) :: y
            if (0.25_RKG < x) then
                y = x + 0.25_RKG
            elseif (FPMIN <= x) then
                y = -0.6931_RKG / log(x)
            else
                y = -0.6931_RKG / log(FPMIN)
            end if
        end function alfa

        pure function exmin1minx(x) result(y)
            ! computes (exp(x) - 1 - x) / (0.5 * x * x).
            real(RKG), intent(in) :: x
            real(RKG) :: y
            real(RKG) :: t, tsq
            if (x == 0) then
                y = 1._RKG
            elseif (0.9_RKG < abs(x)) then
                y = (exp(x) - 1._RKG - x) / (x * x / 2._RKG)
            else
                t = sinh(.5_RKG * x)
                tsq = t**2
                y = (2._RKG * tsq + (2._RKG * t * sqrt(1._RKG + tsq) - x)) / (.5_RKG * x * x)
            end if
        end function exmin1minx

        pure function lngam1(x)
            real(RKG), intent(in) :: x
            real(RKG) :: lngam1
            ! ln(gamma(1+x)), -1<=x<=1
            lngam1 = -logoneplusx(x * (x - 1._RKG) * auxgam(x))
        end function lngam1

        pure function logoneplusx(x) result(y)
            ! x > -1: computes ln(1 + x) with good relative precision when |x| is small.
            real(RKG), intent(in) :: x
            real(RKG) :: y, r, s
            y = log(1._RKG + x)
            if (-0.2928_RKG < x .and. x < 0.4142_RKG) then
                s = y * exmin1(y)
                r = (s - x) / (s + 1._RKG)
                y = y - r * (6._RKG - r) / (6._RKG - 4._RKG * r)
            end if
        end function logoneplusx

        pure function dompart(a, x, qt)
            ! dompart is approx. of  x^a * exp(-x) / gamma(a+1).
            logical(LK), intent(in) :: qt
            real(RKG), intent(in) :: a, x
            real(RKG) :: dompart, logx, c, dp, la, mu, r
            logx = log(x)
            if (a <= 1._RKG) then
                r = -x + a * logx
            else
                if (x == a) then
                    r = 0._RKG
                else
                    la = x / a
                    r = a * (1._RKG - la + log(la))
                end if
                r = r - 0.5_RKG * log(6.2832_RKG * a)
            end if
            if (r < LOG_FPMIN) then ! explow
                dp = 0._RKG
            else
                dp = exp(r)
            end if
            if (qt) then
                dompart = dp
            else
                if (a < 3 .or. x < 0.2_RKG) then
                    dompart = exp(a * logx - x) / gamma(a + 1._RKG)
                else
                    mu = (x - a) / a
                    c = lnec(mu)
                    if (a * c > log(FPMAX)) then
                        dompart = -100._RKG
                    else
                        dompart = exp(a * c) / (sqrt(a * 2 * PI) * gamstar(a))
                    end if
                end if
            end if
        end function dompart

        pure function gamstar(x)
            ! {gamstar(x)=exp(stirling(x)), x>0 or }
            ! {gamma(x)/(exp(-x+(x-0.5)*ln(x))/sqrt(2pi)}
            real(RKG), intent(in) :: x
            real(RKG) :: gamstar
            if (x >= 3._RKG) then
                gamstar = exp(stirling(x))
            elseif (x > 0._RKG) then
                gamstar = gamma(x) / (exp(-x + (x - 0.5_RKG) * log(x)) * SQRT_2PI)
            else
                gamstar = FPMAX
            end if
        end function gamstar

        pure function stirling(x)
            !{stirling series, function corresponding with}
            !{asymptotic series for log(gamma(x))}
            !{that is:  1/(12x)-1/(360x**3)... x>= 3}
            real(RKG), intent(in) :: x
            real(RKG) :: stirling, a(0:17), c(0:6), z
            if (x < FPMIN) then
                stirling = FPMAX
            elseif (x < 1._RKG) then
                stirling = lngam1(x) - (x + 0.5_RKG) * log(x) + x - LOG_SQRT_2PI
            elseif (x < 2._RKG) then
                stirling = lngam1(x - 1) - (x - 0.5_RKG) * log(x) + x - LOG_SQRT_2PI
            elseif (x < 3._RKG) then
                stirling = lngam1(x - 2) - (x - 0.5_RKG) * log(x) + x - LOG_SQRT_2PI + log(x - 1._RKG)
            elseif (x < 12._RKG) then
                a(0) = 1.996379051590076518221_RKG
                a(1) = -0.17971032528832887213e-2_RKG
                a(2) = 0.131292857963846713e-4_RKG
                a(3) = -0.2340875228178749e-6_RKG
                a(4) = 0.72291210671127e-8_RKG
                a(5) = -0.3280997607821e-9_RKG
                a(6) = 0.198750709010e-10_RKG
                a(7) = -0.15092141830e-11_RKG
                a(8) = 0.1375340084e-12_RKG
                a(9) = -0.145728923e-13_RKG
                a(10) = 0.17532367e-14_RKG
                a(11) = -0.2351465e-15_RKG
                a(12) = 0.346551e-16_RKG
                a(13) = -0.55471e-17_RKG
                a(14) = 0.9548e-18_RKG
                a(15) = -0.1748e-18_RKG
                a(16) = 0.332e-19_RKG
                a(17) = -0.58e-20_RKG
                z = 18._RKG / (x * x) - 1._RKG
                stirling = chepolsum(17, z, a) / (12._RKG * x)
            else
                z = 1._RKG / (x * x)
                if (x < 1000._RKG) then
                    c(0) = 0.25721014990011306473e-1_RKG
                    c(1) = 0.82475966166999631057e-1_RKG
                    c(2) = -0.25328157302663562668e-2_RKG
                    c(3) = 0.60992926669463371e-3_RKG
                    c(4) = -0.33543297638406e-3_RKG
                    c(5) = 0.250505279903e-3_RKG
                    c(6) = 0.30865217988013567769_RKG
                    stirling = ((((((c(5) * z + c(4)) * z + c(3)) * z + c(2)) * z + c(1)) * z + c(0)) / (c(6) + z) / x)
                else
                    stirling = (((-z / 1680._RKG + 1._RKG / 1260._RKG) * z - 1._RKG / 360._RKG) * z + 1._RKG / 12._RKG) / x
                end if
            end if
        end function stirling

        pure function pqasymp(a, x, dp, p)
            logical(LK), intent(in) :: p
            real(RKG), intent(in) :: a, x, dp
            real(RKG) :: pqasymp, y, mu, eta, u, v
            integer(IK) :: s
            if (dp == 0._RKG) then
                if (p) then
                    pqasymp = 0._RKG
                else
                    pqasymp = 1._RKG
                end if
            else
                if (p) then
                    s = -1_IK
                else
                    s = 1_IK
                end if
                mu = (x - a) / a
                y = -lnec(mu)
                if (y < 0._RKG) then
                    eta = 0._RKG
                else
                    eta = sqrt(2._RKG * y)
                end if
                y = y * a
                v = sqrt(abs(y))
                if (mu < 0._RKG) then
                    eta = -eta
                    v = -v
                end if
                ! warning: `errorfunction` custom function replaced with Fortran intrinsic function `erfc()`.
                u = 0.5_RKG * erfc(s * v) ! errorfunction(s * v, .true., .false.)
                v = s * exp(-y) * saeta(a, eta) / sqrt(2._RKG * PI * a)
                pqasymp = u + v
            end if
        end function pqasymp

        pure function saeta(a, eta)
            real(RKG), intent(in) :: a, eta
            real(RKG) :: saeta, y, s, t, fm(0 : 26), bm(0 : 26)
            integer(IK) :: m
            fm(0) = 1.0_RKG
            fm(1) = -1.0_RKG/3.0_RKG
            fm(2) = 1.0_RKG/12.0_RKG
            fm(3) = -2.0_RKG/135.0_RKG
            fm(4) = 1.0_RKG/864.0_RKG
            fm(5) = 1.0_RKG/ 2835.0_RKG
            fm(6) = -139.0_RKG/777600.0_RKG
            fm(7) = 1.0_RKG/25515.0_RKG
            fm(8) = -571.0_RKG/261273600.0_RKG
            fm(9) = -281.0_RKG/151559100.0_RKG
            fm(10) = 8.29671134095308601e-7_RKG
            fm(11) = -1.76659527368260793e-7_RKG
            fm(12) = 6.70785354340149857e-9_RKG
            fm(13) = 1.02618097842403080e-8_RKG
            fm(14) = -4.38203601845335319e-9_RKG
            fm(15) = 9.14769958223679023e-10_RKG
            fm(16) = -2.55141939949462497e-11_RKG
            fm(17) = -5.83077213255042507e-11_RKG
            fm(18) = 2.43619480206674162e-11_RKG
            fm(19) = -5.02766928011417559e-12_RKG
            fm(20) = 1.10043920319561347e-13_RKG
            fm(21) = 3.37176326240098538e-13_RKG
            fm(22) = -1.39238872241816207e-13_RKG
            fm(23) = 2.85348938070474432e-14_RKG
            fm(24) = -5.13911183424257258e-16_RKG
            fm(25) = -1.97522882943494428e-15_RKG
            fm(26) =  8.09952115670456133e-16_RKG
            bm(25) = fm(26)
            bm(24) = fm(25)
            do m = 24_IK, 1_IK, -1_IK
                bm(m - 1_IK) = fm(m) + (m + 1_IK) * bm(m + 1_IK) / a
            end do
            s = bm(0)
            t = s
            y = eta
            m = 1_IK
            do while (abs(t / s) > EPS .and. m < 25_IK)
                t = bm(m) * y
                s = s + t
                m = m + 1_IK
                y = y * eta
            end do
            saeta = s / (1._RKG + bm(1) / a)
        end function saeta

        pure function qfraction(kappa, x, dp) result(gammaIncUpp)
            real(RKG), intent(in) :: kappa, x, dp
            real(RKG) :: g, gammaIncLow, gammaIncUpp, r, s, t, tau, ro
            if (dp == 0._RKG) then
                gammaIncUpp = 0._RKG
            else
                gammaIncLow = 0
                gammaIncUpp = (x - 1._RKG - kappa) * (x + 1._RKG - kappa)
                r = 4 * (x + 1._RKG - kappa)
                s = 1._RKG - kappa
                ro = 0._RKG
                t = 1.0_RKG
                g = 1.0_RKG
                do while (EPS <= abs(t / g))
                    gammaIncLow = gammaIncLow + s
                    gammaIncUpp = gammaIncUpp + r
                    r = r + 8._RKG
                    s = s + 2._RKG
                    tau = gammaIncLow * (1._RKG + ro)
                    ro = tau / (gammaIncUpp - tau)
                    t = ro * t
                    g = g + t
                end do
                gammaIncUpp = (kappa / (x + 1._RKG - kappa)) * g * dp
            end if
        end function qfraction

        pure function qtaylor(kappa, x, dp) result(gammaIncUpp)
            real(RKG), intent(in) :: kappa, x, dp
            real(RKG) :: logx, gammaIncLow, gammaIncUpp, r, s, t, u, v
            logx = log(x)
            if (dp == 0._RKG) then
                gammaIncUpp = 0._RKG
            else
                r = kappa * logx
                gammaIncUpp = r * exmin1(r) ! {gammaIncUpp = x^kappa - 1 }
                s = kappa * (1._RKG - kappa) * auxgam(kappa) ! {s = 1-1/gamma(1+kappa) }
                gammaIncUpp = (1 - s) * gammaIncUpp
                u = s - gammaIncUpp ! {u = 1 - x^kappa/gamma(1+kappa)}
                gammaIncLow = kappa * x
                gammaIncUpp = kappa + 1._RKG
                r = kappa + 3._RKG
                t = 1._RKG
                v = 1._RKG
                do while (abs(t / v) > EPS)
                    gammaIncLow = gammaIncLow + x
                    gammaIncUpp = gammaIncUpp + r
                    r = r + 2._RKG
                    t = -gammaIncLow * t / gammaIncUpp
                    v = v + t
                end do
                v = kappa * (1._RKG - s) * exp((kappa + 1._RKG) * logx) * v / (kappa + 1._RKG)
                gammaIncUpp = u + v
            end if
        end function qtaylor

        pure function exmin1(x) result(y)
            ! computes (exp(x) - 1) / x.
            real(RKG), intent(in) :: x
            real(RKG) :: t, y
            if (x == 0._RKG) then
                y = 1._RKG
            elseif (x < -0.69_RKG .or. 0.4_RKG < x) then
                y = (exp(x) - 1._RKG) / x
            else
                t = x / 2._RKG
                y = exp(t) * sinh(t) / t
            end if
        end function exmin1

        pure recursive function auxgam(x) result(auxgamm)
            ! function g in 1/gamma(x+1)=1+x*(x-1)*g(x), -1 <= x <= 1.
            real(RKG), intent(in) :: x
            real(RKG) :: auxgamm, t, dr(0:17)
            if (x < 0._RKG) then
                auxgamm = -(1._RKG + (1._RKG + x) * (1._RKG + x) * auxgam(1._RKG + x)) / (1._RKG - x)
            else
                dr(0)   = -1.013609258009865776949_RKG
                dr(1)   = 0.784903531024782283535e-1_RKG
                dr(2)   = 0.67588668743258315530e-2_RKG
                dr(3)   = -0.12790434869623468120e-2_RKG
                dr(4)   = 0.462939838642739585e-4_RKG
                dr(5)   = 0.43381681744740352e-5_RKG
                dr(6)   = -0.5326872422618006e-6_RKG
                dr(7)   = 0.172233457410539e-7_RKG
                dr(8)   = 0.8300542107118e-9_RKG
                dr(9)   = -0.10553994239968e-9_RKG
                dr(10)  = 0.39415842851e-11_RKG
                dr(11)  = 0.362068537e-13_RKG
                dr(12)  = -0.107440229e-13_RKG
                dr(13)  = 0.5000413e-15_RKG
                dr(14)  = -0.62452e-17_RKG
                dr(15)  = -0.5185e-18_RKG
                dr(16)  = 0.347e-19_RKG
                dr(17)  = -0.9e-21_RKG
                t = 2 * x - 1._RKG
                auxgamm = chepolsum(17_IK, t, dr)
            end if
        end function auxgam

        pure function chepolsum(n, x, a)
            real(RKG), intent(in) :: x, a(0:n)
            integer(IK), intent(in) :: n
            real(RKG) :: chepolsum
            real(RKG) :: h, r, s, tx
            integer(IK) :: k
            ! a[0]/2+a[1]t1(x)+...a[n]tn(x) series of chebychev polynomials.
            if (n == 0_IK) then
                chepolsum = a(0) / 2._RKG
            elseif (n == 1_IK) then
                chepolsum = a(0) / 2._RKG + a(1) * x
            else
                tx = x + x
                r = a(n)
                h = a(n - 1_IK) + r * tx
                do k = n - 2_IK, 1_IK, -1_IK
                    s = r
                    r = h
                    h = a(k) + r * tx - s
                end do
                chepolsum = a(0) / 2._RKG - r + h * x
            end if
        end function chepolsum

        pure function ptaylor(kappa, x, dp) result(gammaIncLow)
            real(RKG), intent(in) :: kappa, x, dp
            real(RKG) :: gammaIncLow, c, r
            if (dp == 0._RKG) then
                gammaIncLow = 0._RKG
            else
                gammaIncLow = 1._RKG
                c = 1._RKG
                r = kappa
                do while ((c / gammaIncLow) > EPS)
                    r = r + 1._RKG
                    c = x * c / r
                    gammaIncLow = gammaIncLow + c
                end do
                gammaIncLow = gammaIncLow * dp
            end if
        end function ptaylor

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif