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
!>  [getGammaIncUpp](@ref pm_mathGamma::getGammaIncUpp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: LK
        use pm_val2str, only: getStr

        real(RK), parameter :: EPS = epsilon(0._RK) * 100
        real(RK), allocatable :: Point(:), Kappa(:), GammaIncUpp(:), GammaIncUpp_ref(:), diff(:), Tol(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [real(RK)::]
        Point = [real(RK)::]
        Kappa = [real(RK)::]
        GammaIncUpp_ref = [real(RK)::]

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `GammaIncUpp` with empty input `Point` and `Kappa`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `GammaIncUpp` with empty input `Point` and `Kappa`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `GammaIncUpp` with empty input `Point` and `Kappa` and `Tol`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [0._RK]
        Kappa = [1.0_RK]
        GammaIncUpp_ref = [1._RK]

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [0.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [0.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [0.0_RK], shape = [1._RK], tol = EPS`.")

        block
            use pm_kind, only: LK
            logical(LK) :: Failed(size(Point))
            call setGammaIncUppContFrac(GammaIncUpp, Point, Kappa, failed, Tol)
            call report()
            call test%assert(assertion, SK_"getGammaIncUppContFrac() must yield a correct `GammaIncUpp` with `Point = [0.0_RK], shape = [1._RK], tol = EPS`.")
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [1.0_RK]
        Kappa = [1.0_RK]
        GammaIncUpp_ref = 1._RK - [0.632120558828557678404476229838539224_RK]

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [3.0_RK]
        Kappa = [5.0_RK]
        GammaIncUpp_ref = 1._RK - [0.184736755476227933713267943730238524_RK]

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [7.0_RK]
        Kappa = [0.5_RK]
        GammaIncUpp_ref = 1._RK - [0.999817189367018164968240003449102634_RK]

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK]`.")

        GammaIncUpp = getGammaIncUpp(Point, Kappa, Tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tol = [EPS]
        Point = [7.0_RK]
        Kappa = [0.5_RK]
        GammaIncUpp_ref = 1._RK - [0.999817189367018164968240003449102634_RK]

        GammaIncUpp = [getGammaIncUpp(Point(1), Kappa(1))]
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = 1.0_RK, shape = 1._RK`.")

        GammaIncUpp = getGammaIncUpp(Point(1), Kappa(1))
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = 1.0_RK, shape = 1._RK`.")

        GammaIncUpp = getGammaIncUpp(Point(1), Kappa(1), Tol(1))
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `GammaIncUpp` with `Point = 1.0_RK, shape = 1._RK, tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(GammaIncUpp - GammaIncUpp_ref)
            assertion = assertion .and. all(abs(diff) <= Tol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Point          ", Point
                write(test%disp%unit,"(*(g0,:,', '))") "GammaIncUpp    ", GammaIncUpp
                write(test%disp%unit,"(*(g0,:,', '))") "GammaIncUpp_ref", GammaIncUpp_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
