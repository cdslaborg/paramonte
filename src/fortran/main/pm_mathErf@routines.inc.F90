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
!>  This include file contains the implementation of procedures in [pm_mathErf](@ref pm_mathErf).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%
#if     getErfInv_ENABLED
        !%%%%%%%%%%%%%%%%

        if (present(abserr)) then
            call setErfInv(erfinv, x, abserr)
        else
            call setErfInv(erfinv, x, 100 * epsilon(0._RKG))
        end if

        !%%%%%%%%%%%%%%%%
#elif   setErfInv_ENABLED
        !%%%%%%%%%%%%%%%%

        real(RKG) :: temp, absx
        real(RKG), parameter :: TWO_OVER_SQRTPI = 2._RKG / sqrt(acos(-1._RKG))
        ! Chebyshev declarations.
        integer(IK) , parameter :: MINDEL = 0, NDELTA = 37, NLAMDA = 26, NMU = 25, NXI = 38
        integer(IK) , parameter :: MAXDEL = MINDEL + NDELTA ! 37
        integer(IK) , parameter :: MINLAM = MAXDEL + 1 ! 38
        integer(IK) , parameter :: MAXLAM = MINLAM + NLAMDA ! 64
        integer(IK) , parameter :: MINMU = MAXLAM + 1 ! 65
        integer(IK) , parameter :: MAXMU = MINMU + NMU ! 90
        integer(IK) , parameter :: MINXI = MAXMU + 1 ! 91
        integer(IK) , parameter :: MAXXI = MINXI + NXI ! 129
        ! The vector `chebyshScaler` is used to scale the argument of the Chebyshev polynomial expansion in the range 1.D-300 < 1 - X < 0.2.
        real(RKG), parameter :: chebyshScaler(6) =  [ -1.548813042373261659512742_RKG &
                                                    , +2.565490123147816151928163_RKG &
                                                    , -.5594576313298323225436913_RKG &
                                                    , +2.287915716263357638965891_RKG &
                                                    , -9.199992358830151031278420_RKG &
                                                    , +2.794990820124599493768426_RKG ]
        ! The coefficients of polynomial expansions, stored in the order DELTA(0..37), LAMDA(0..26), MU(0..25), XI(0..38), where
        !   DELTA are coefficients of the Chebyshev polynomial expansion of erfinv(X) for 0.9975 < X <= 1-5.0D-16.
        !   LAMDA are coefficients of the Chebyshev polynomial expansion of erfinv(X) for 0.8 < X <= 0.9975.
        !   MU    are coefficients of the Chebyshev polynomial expansion of erfinv(X) for 1 - 5.0D-16 < X <= 1 - 1.D-300.
        !   XI    are coefficients of the Chebyshev polynomial expansion of erfinv(X) in the range 0.0 <= X <= 0.8.
        real(RKG), parameter :: coef(0:MAXXI) = [ +0.9566797090204925274526373_RKG &
                                                , -0.0231070043090649036999908_RKG &
                                                , -0.0043742360975084077333218_RKG &
                                                , -0.0005765034226511854809364_RKG &
                                                , -0.0000109610223070923931242_RKG &
                                                , +0.0000251085470246442787982_RKG &
                                                , +0.0000105623360679477511955_RKG &
                                                , +0.0000027544123300306391503_RKG &
                                                , +0.0000004324844983283380689_RKG &
                                                , -0.0000000205303366552086916_RKG &
                                                , -0.0000000438915366654316784_RKG &
                                                , -0.0000000176840095080881795_RKG &
                                                , -0.0000000039912890280463420_RKG &
                                                , -0.0000000001869324124559212_RKG &
                                                , +0.0000000002729227396746077_RKG &
                                                , +0.0000000001328172131565497_RKG &
                                                , +0.0000000000318342484482286_RKG &
                                                , +0.0000000000016700607751926_RKG &
                                                , -0.0000000000020364649611537_RKG &
                                                , -0.0000000000009648468127965_RKG &
                                                , -0.0000000000002195672778128_RKG &
                                                , -0.0000000000000095689813014_RKG &
                                                , +0.0000000000000137032572230_RKG &
                                                , +0.0000000000000062538505417_RKG &
                                                , +0.0000000000000014584615266_RKG &
                                                , +0.0000000000000001078123993_RKG &
                                                , -0.0000000000000000709229988_RKG &
                                                , -0.0000000000000000391411775_RKG &
                                                , -0.0000000000000000111659209_RKG &
                                                , -0.0000000000000000015770366_RKG &
                                                , +0.0000000000000000002853149_RKG &
                                                , +0.0000000000000000002716662_RKG &
                                                , +0.0000000000000000000957770_RKG &
                                                , +0.0000000000000000000176835_RKG &
                                                , -0.0000000000000000000009828_RKG &
                                                , -0.0000000000000000000020464_RKG &
                                                , -0.0000000000000000000008020_RKG &
                                                , -0.0000000000000000000001650_RKG &
                                                , +0.9121588034175537733059200_RKG &
                                                , -0.0162662818676636958546661_RKG &
                                                , +0.0004335564729494453650589_RKG &
                                                , +0.0002144385700744592065205_RKG &
                                                , +0.0000026257510757648130176_RKG &
                                                , -0.0000030210910501037969912_RKG &
                                                , -0.0000000124060618367572157_RKG &
                                                , +0.0000000624066092999917380_RKG &
                                                , -0.0000000005401247900957858_RKG &
                                                , -0.0000000014232078975315910_RKG &
                                                , +0.0000000000343840281955305_RKG &
                                                , +0.0000000000335848703900138_RKG &
                                                , -0.0000000000014584288516512_RKG &
                                                , -0.0000000000008102174258833_RKG &
                                                , +0.0000000000000525324085874_RKG &
                                                , +0.0000000000000197115408612_RKG &
                                                , -0.0000000000000017494333828_RKG &
                                                , -0.0000000000000004800596619_RKG &
                                                , +0.0000000000000000557302987_RKG &
                                                , +0.0000000000000000116326054_RKG &
                                                , -0.0000000000000000017262489_RKG &
                                                , -0.0000000000000000002784973_RKG &
                                                , +0.0000000000000000000524481_RKG &
                                                , +0.0000000000000000000065270_RKG &
                                                , -0.0000000000000000000015707_RKG &
                                                , -0.0000000000000000000001475_RKG &
                                                , +0.0000000000000000000000450_RKG &
                                                , +0.9885750640661893136460358_RKG &
                                                , +0.0108577051845994776160281_RKG &
                                                , -0.0017511651027627952494825_RKG &
                                                , +0.0000211969932065633437984_RKG &
                                                , +0.0000156648714042435087911_RKG &
                                                , -0.0000005190416869103124261_RKG &
                                                , -0.0000000371357897426717780_RKG &
                                                , +0.0000000012174308662357429_RKG &
                                                , -0.0000000001768115526613442_RKG &
                                                , -0.0000000000119372182556161_RKG &
                                                , +0.0000000000003802505358299_RKG &
                                                , -0.0000000000000660188322362_RKG &
                                                , -0.0000000000000087917055170_RKG &
                                                , -0.0000000000000003506869329_RKG &
                                                , -0.0000000000000000697221497_RKG &
                                                , -0.0000000000000000109567941_RKG &
                                                , -0.0000000000000000011536390_RKG &
                                                , -0.0000000000000000001326235_RKG &
                                                , -0.0000000000000000000263938_RKG &
                                                , +0.0000000000000000000005341_RKG &
                                                , -0.0000000000000000000022610_RKG &
                                                , +0.0000000000000000000009552_RKG &
                                                , -0.0000000000000000000005250_RKG &
                                                , +0.0000000000000000000002487_RKG &
                                                , -0.0000000000000000000001134_RKG &
                                                , +0.0000000000000000000000420_RKG &
                                                , +0.9928853766189408231495800_RKG &
                                                , +0.1204675161431044864647846_RKG &
                                                , +0.0160781993420999447257039_RKG &
                                                , +0.0026867044371623158279591_RKG &
                                                , +0.0004996347302357262947170_RKG &
                                                , +0.0000988982185991204409911_RKG &
                                                , +0.0000203918127639944337340_RKG &
                                                , +0.0000043272716177354218758_RKG &
                                                , +0.0000009380814128593406758_RKG &
                                                , +0.0000002067347208683427411_RKG &
                                                , +0.0000000461596991054300078_RKG &
                                                , +0.0000000104166797027146217_RKG &
                                                , +0.0000000023715009995921222_RKG &
                                                , +0.0000000005439284068471390_RKG &
                                                , +0.0000000001255489864097987_RKG &
                                                , +0.0000000000291381803663201_RKG &
                                                , +0.0000000000067949421808797_RKG &
                                                , +0.0000000000015912343331469_RKG &
                                                , +0.0000000000003740250585245_RKG &
                                                , +0.0000000000000882087762421_RKG &
                                                , +0.0000000000000208650897725_RKG &
                                                , +0.0000000000000049488041039_RKG &
                                                , +0.0000000000000011766394740_RKG &
                                                , +0.0000000000000002803855725_RKG &
                                                , +0.0000000000000000669506638_RKG &
                                                , +0.0000000000000000160165495_RKG &
                                                , +0.0000000000000000038382583_RKG &
                                                , +0.0000000000000000009212851_RKG &
                                                , +0.0000000000000000002214615_RKG &
                                                , +0.0000000000000000000533091_RKG &
                                                , +0.0000000000000000000128488_RKG &
                                                , +0.0000000000000000000031006_RKG &
                                                , +0.0000000000000000000007491_RKG &
                                                , +0.0000000000000000000001812_RKG &
                                                , +0.0000000000000000000000439_RKG &
                                                , +0.0000000000000000000000106_RKG &
                                                , +0.0000000000000000000000026_RKG &
                                                , +0.0000000000000000000000006_RKG &
                                                , +0.0000000000000000000000002_RKG &
                                                ]
        real(RKG), parameter :: HALF_EPS = 0.5_RKG * epsilon(0._RKG)
        real(RKG) :: carg, cargTwice, w1, w2, w3
        ! COEF_LIM_FULL is an array containing MINXI, MAXXI, MINLAM, MAXLAM, MINDEL, MAXDEL, MINMU, MAXMU in locations -1..6 where,
        !       MAX...  is the position in `coef` of the last coefficient of a Chebyshev polynomial expansion.
        !       MIN...  is the position in `coef` of the first coefficient of a Chebyshev polynomial expansion.
        integer(IK) , parameter :: COEF_LIM_FULL(-1:6) = [MINXI, MAXXI, MINLAM, MAXLAM, MINDEL, MAXDEL, MINMU, MAXMU]
        ! Decide which approximation to use, and calculate the argument of the Chebyshev polynomial expansion.
        ! Compute the minimum index of a coefficient in the Chebyshev polynomial expansion to be used.
        ! That is, include only the coefficients whose precision is relevant to the requested output precision.
        ! \bug ifort 2021.8 : Cannot digest implied do loop below. Gfortran is fine.
        ! \remedy For now, the loop was unrolled.
       !integer(IK) , parameter :: COEF_LIM(-1:6) = [(COEF_LIM_FULL(i - 1), COEF_LIM_FULL(i) - count(abs(coef(COEF_LIM_FULL(i - 1) : COEF_LIM_FULL(i))) < HALF_EPS), i = 0, 6, 2)]
        integer(IK) , parameter :: COEF_LIM(-1:6) = [ COEF_LIM_FULL(-1), COEF_LIM_FULL(0) - count(abs(coef(COEF_LIM_FULL(-1) : COEF_LIM_FULL(0))) < HALF_EPS) &
                                                    , COEF_LIM_FULL(+1), COEF_LIM_FULL(2) - count(abs(coef(COEF_LIM_FULL(+1) : COEF_LIM_FULL(2))) < HALF_EPS) &
                                                    , COEF_LIM_FULL(+3), COEF_LIM_FULL(4) - count(abs(coef(COEF_LIM_FULL(+3) : COEF_LIM_FULL(4))) < HALF_EPS) &
                                                    , COEF_LIM_FULL(+5), COEF_LIM_FULL(6) - count(abs(coef(COEF_LIM_FULL(+5) : COEF_LIM_FULL(6))) < HALF_EPS) &
                                                    ]
        integer(IK) :: i, j, imin
        if(-1._RKG < x .and. x < 1._RKG) then
            if(1.e-7_RKG < abserr) then
                temp = -log((1.0 - x)*(1.0 + x));
                if (temp < 5._RKG) then
                    temp = temp - 2.5_RKG
                    erfinv = +2.81022636e-08_RKG
                    erfinv = +3.43273939e-07_RKG + erfinv * temp
                    erfinv = -3.5233877e-06_RKG + erfinv * temp
                    erfinv = -4.39150654e-06_RKG + erfinv * temp
                    erfinv = +0.00021858087_RKG + erfinv * temp
                    erfinv = -0.00125372503_RKG + erfinv * temp
                    erfinv = -0.00417768164_RKG + erfinv * temp
                    erfinv = +0.246640727_RKG + erfinv * temp
                    erfinv = +1.50140941_RKG + erfinv * temp
                else
                    temp = sqrt(temp) - 3._RKG
                    erfinv = -0.000200214257_RKG
                    erfinv = +0.000100950558_RKG + erfinv * temp
                    erfinv = +0.00134934322_RKG + erfinv * temp
                    erfinv = -0.00367342844_RKG + erfinv * temp
                    erfinv = +0.00573950773_RKG + erfinv * temp
                    erfinv = -0.0076224613_RKG + erfinv * temp
                    erfinv = +0.00943887047_RKG + erfinv * temp
                    erfinv = +1.00167406_RKG + erfinv * temp
                    erfinv = +2.83297682_RKG + erfinv * temp
                end if
                erfinv = erfinv * x
            elseif (2.e-24_RKG < abserr) then
                !	This algorithm is based on the approximate formulae of Mathematics of Computation 22, (1968) PP144-158, which he claims to be accurate up to 23 decimal points.
                !	Amir Shahmoradi, Nov 10, 2009, 8:53 PM, Michigan
                !	Amir Shahmoradi, Wednesday 5:30 AM, XXX XX, 2012, Institute for Fusion Studies, The University of Texas Austin<br>
                !	HAPPY BIRTHDAY AMIR! (Self-congrats syndrome caused by lack of sleep)
                if (x == 0._RKG) then
                    erfinv = 0._RKG
                    return
                end if
                absx = abs(x) ! The argument of the Chebyshev polynomial expansion.
                if (absx <= 0.8_RKG) then
                    carg = 3.125_RKG * absx * absx - 1._RKG
                    j = -1
                else
                    if (absx <= 0.9975_RKG) then
                        j = 1
                    else
                        j = 3
                    end if
                    absx = sqrt(-log((1._RKG - absx) * (1._RKG + absx)))
                    carg = chebyshScaler(j) * absx + chebyshScaler(j + 1)
                end if
                ! Evaluate the Chebyshev polynomial expansion.
                cargTwice = carg + carg
                ! W1..W3 are the adjacent elements of the recurrence used to evaluate the Chebyshev polynomial expansion.
                w1 = 0._RKG
                w2 = 0._RKG
                imin = COEF_LIM(j)
                i = COEF_LIM(j + 1)
                do
                    w3 = w2
                    w2 = w1
                    w1 = (cargTwice * w2 - w3) + coef(i)
                    i = i - 1
                    if (imin < i) cycle
                    exit
                end do
                erfinv = sign(absx * ((carg * w1 - w2) + coef(imin)), x)
            else
                ! Use Halley approximation. 
                ! Verified accuracy up to at least 2 * 10**-26.
                absx = 1._RKG - abs(x)
                temp = sqrt(-2._RKG * log(0.5_RKG * absx))
                erfinv = -0.70711_RKG * ((2.30753_RKG + temp * 0.27061_RKG) / (1._RKG + temp * (0.99229_RKG + temp * 0.04481_RKG)) - temp)
                do i = 1, 2
                    temp = erfc(erfinv) - absx
                    erfinv = erfinv + temp / (TWO_OVER_SQRTPI * exp(-erfinv**2) - erfinv * temp)
                end do
                if (x < 0._RKG) erfinv = -erfinv
            end if
        else
            CHECK_ASSERTION(__LINE__, -1._RKG < x .and. x < 1._RKG, SK_"@setErfInv(): The condition `-1. < x .and. x < 1.` must hold. x = "//getStr(x)) ! fpp
            erfinv = sign(huge(x), x)
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  strecok_ENABLED
#undef  halley_ENABLED