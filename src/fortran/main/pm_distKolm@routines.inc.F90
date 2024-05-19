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
!>  This include file contains the implementation of procedures in [pm_distKolm](@ref pm_distKolm).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getKolmPDF_ENABLED || setKolmPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: niter, niterSq
        real(TKG) :: expsq, temp, tempold, delta, xinv, exponent
        real(TKG), parameter :: xbreak = 1._TKG ! this seems to be the performance tipping point between the two series convergence in real64 and real32 precision.
        real(TKG), parameter :: SQRT2 = sqrt(2._TKG)
        real(TKG), parameter :: TWOSQRT2PI = 2 * SQRT2 * sqrt(acos(-1._TKG))
        real(TKG), parameter :: PI2SQRT8 = acos(-1._TKG) / sqrt(8._TKG)
        real(TKG), parameter :: TOL = 10 * epsilon(0._TKG)
        niter = 1_IK
        ! \devnote
        ! Based on experimentations, it takes roughly 2-4 term of the series in each branch to compute the PDF in [0, 3].
        if (xbreak < x) then
            exponent = SQRT2 * x
            expsq = -exponent**2
            pdf = exp(expsq) ! first term in the series.
            if (pdf + TOL < 1._TKG) then
                tempold = pdf
                delta = 1._TKG
                do
                    niter = niter + 1_IK
                    niterSq = niter * niter
                    temp = niterSq * exp(expsq * niterSq)
                    delta = sign(temp, -delta)
                    pdf = pdf + delta
                    if (tempold - temp < TOL) exit
                    tempold = temp
                end do
                pdf = 8 * x * pdf
            else
                pdf = 0._TKG
            end if
            !print *, (niter + 1_IK) / 2_IK
        elseif (0._TKG < x) then ! warning: `x < 1` must hold.
            xinv = 1._TKG / x
            exponent = PI2SQRT8 * xinv
            expsq = exponent**2
            pdf = (expsq - .5_TKG) * exp(-expsq)
            do
                niter = niter + 2_IK
                expsq = (exponent * niter)**2
                delta = (expsq - .5_TKG) * exp(-expsq)
                pdf = pdf + delta
                if (delta < TOL) exit
            end do
            pdf = pdf * TWOSQRT2PI * xinv**2
            !print *, (niter + 1_IK) / 2_IK
        else
            CHECK_ASSERTION(__LINE__, 0._TKG <= x, SK_"@setKolmPDF(): The condition `0 < x` must hold. x = "//getStr(x))
            pdf = 0._TKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getKolmCDF_ENABLED || setKolmCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: niter
        real(TKG) :: temp, tempold, delta, xinv, exponent
        real(TKG), parameter :: xbreak = 1._TKG ! this seems to be the performance tipping point between the two series convergence in real64 and real32 precision.
        real(TKG), parameter :: SQRT2 = sqrt(2._TKG)
        real(TKG), parameter :: SQRT2PI = SQRT2 * sqrt(acos(-1._TKG))
        real(TKG), parameter :: PI2SQRT8 = acos(-1._TKG) / sqrt(8._TKG)
        real(TKG), parameter :: TOL = 10 * epsilon(0._TKG)
        niter = 1_IK
        ! \devnote
        ! Based on experimentations, it takes roughly two term of the series in each branch to compute the CDF.
        if (xbreak < x) then
            exponent = SQRT2 * x
            cdf = exp(-exponent**2) ! first term in the series.
            if (cdf + TOL < 1._TKG) then
                tempold = cdf
                delta = 1._TKG
                do
                    niter = niter + 1_IK
                    temp = exp(-(exponent * niter)**2)
                    delta = sign(temp, -delta)
                    cdf = cdf + delta
                    if (tempold - temp < TOL) exit
                    tempold = temp
                end do
                cdf = 1._TKG - 2 * cdf
            else
                cdf = 0._TKG
            end if
            !print *, (niter + 1_IK) / 2_IK
        elseif (0._TKG < x) then
            xinv = 1._TKG / x
            exponent = PI2SQRT8 * xinv
            cdf = exp(-exponent**2)
            do
                niter = niter + 2_IK
                delta = exp(-(exponent * niter)**2)
                cdf = cdf + delta
                if (delta < TOL) exit
            end do
            cdf = cdf * SQRT2PI * xinv
            !print *, (niter + 1_IK) / 2_IK
        else
            CHECK_ASSERTION(__LINE__, 0._TKG <= x, SK_"@setKolmPDF(): The condition `0 < x` must hold. x = "//getStr(x))
            cdf = 0._TKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getKolmQuan_ENABLED || setKolmQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: cdfc
        real(TKG) :: logQuan, quanOld, func, ff, u, t
        real(TKG), parameter :: TOL = 100 * epsilon(0._TKG)
        real(TKG), parameter :: HALFPI = acos(-1._TKG) / 2._TKG
        real(TKG), parameter :: PIOVER8 = HALFPI / 4._TKG
        CHECK_ASSERTION(__LINE__, 0._TKG <= cdf .and. cdf < 1._TKG, SK_"@setKolmQuan(): The condition `0 <= cdf .and. cdf < 1` must hold. cdf = "//getStr(cdf))
        cdfc = 1._TKG - cdf
        if (cdfc < 0.3_TKG) then
            quan = 0.03_TKG
            do
                quanOld = quan
                quan = 0.5_TKG * cdfc + quan**4 - quan**9
                if (quan > 0.06_TKG) quan = quan + quan**16 - quan**25
                if (abs((quanOld - quan) / quan) < TOL) exit
            end do
            quan = sqrt(-0.5_TKG * log(quan))
        elseif (cdfc < 1._TKG) then
            func = -PIOVER8 * (1._TKG - cdfc)**2
            quan = getInvXlogX(func)
            do
                logQuan = log(quan)
                ff = func / (1._TKG + quan**4 + quan**12)**2
                u = (quan * logQuan - ff) / (1._TKG + logQuan)
                t = u / max(0.5_TKG, 1._TKG - 0.5_TKG * u / (quan * (1._TKG + logQuan)))
                quan = quan - t
                if (abs(t / quan) < TOL) exit
            end do
            quan = HALFPI / sqrt(-log(quan))
        else
            quan = 0._TKG
        end if

    contains

        PURE function getInvXlogX(y) result(invXlogX)
            real(TKG), parameter :: TOL = 10 * epsilon(0._TKG)
            real(TKG), parameter :: SQRTOL = sqrt(TOL)
            real(TKG), parameter :: STSTOL = sqrt(SQRTOL)
            real(TKG), parameter :: invNeper = exp(-1._TKG) ! 0.367879441171442322
            real(TKG), intent(in) :: y
            real(TKG) :: invXlogX
            real(TKG) :: t, to
            CHECK_ASSERTION(__LINE__, -invNeper < y .and. y < 0._TKG, SK_"@setKolmQuan(): The condition `-exp(-1) < y .and. y < 0` must hold. y = "//getStr(y))
            to = 0._TKG
            if (y < -0.2) then
                invXlogX = log(invNeper - sqrt(2 * invNeper * (y + invNeper)))
            else
                invXlogX = -10._TKG
            end if
            do
                t = (log(y / invXlogX) - invXlogX) * (invXlogX / (1._TKG + invXlogX))
                invXlogX = invXlogX + t
                if (t < SQRTOL .and. abs(t + to) < STSTOL * abs(t)) exit
                if (abs(t / invXlogX) <= TOL) exit
                to = t
            end do
            invXlogX = exp(invXlogX)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getKolmRand_ENABLED || setKolmRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._TKG <= unif .and. unif < 1._TKG, SK_"@setKolmQuan(): The condition `0 <= unif .and. unif < 1` must hold. unif = "//getStr(unif))
        call setKolmQuan(rand, unif)
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif