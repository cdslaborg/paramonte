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
!>  This include file contains the implementations of the tests of procedures with generic interfaces
!>  [getGammaIncLow](@ref pm_mathGamma::getGammaIncLow).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: LK
        use pm_val2str, only: getStr

        real(RK), parameter :: EPS = epsilon(0._RK) * 100
        real(RK), allocatable :: Point(:), Kappa(:), GammaIncLow(:), GammaIncLow_ref(:), diff(:), Tol(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [real(RK)::]
        Point = [real(RK)::]
        Kappa = [real(RK)::]
        GammaIncLow_ref = [real(RK)::]

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield an empty `GammaIncLow` with empty input `Point` and `Kappa`.")

        GammaIncLow = getGammaIncLow(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield an empty `GammaIncLow` with empty input `Point` and `Kappa` and `Tol`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [0._RK]
        Kappa = [1.0_RK]
        GammaIncLow_ref = [0._RK]

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [0.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [0.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [0.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [1.0_RK]
        Kappa = [1.0_RK]
        GammaIncLow_ref = [0.632120558828557678404476229838539224_RK]

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [3.0_RK]
        Kappa = [5.0_RK]
        GammaIncLow_ref = [0.184736755476227933713267943730238524_RK]

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [7.0_RK]
        Kappa = [0.5_RK]
        GammaIncLow_ref = [0.999817189367018164968240003449102634_RK]

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncLow = getGammaIncLow(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [7.0_RK]
        Kappa = [0.5_RK]
        GammaIncLow_ref = [0.999817189367018164968240003449102634_RK]

        GammaIncLow = [getGammaIncLow(Point(1), Kappa(1))]
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = 1.0_RK, shape = 1._RK`.")

        GammaIncLow = getGammaIncLow(Point(1), Kappa(1), Tol(1))
        call report()
        call test%assert(assertion, SK_"getGammaIncLow() must yield a correct `GammaIncLow` with `Point = 1.0_RK, shape = 1._RK, tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(GammaIncLow - GammaIncLow_ref)
            assertion = assertion .and. all(abs(diff) <= Tol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Point          ", Point
                write(test%disp%unit,"(*(g0,:,', '))") "GammaIncLow    ", GammaIncLow
                write(test%disp%unit,"(*(g0,:,', '))") "GammaIncLow_ref", GammaIncLow_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
