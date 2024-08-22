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
        real(RK), allocatable :: point(:), kappa(:), gammaIncUpp(:), GammaIncUpp_ref(:), diff(:), tol(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [real(RK)::]
        point = [real(RK)::]
        kappa = [real(RK)::]
        GammaIncUpp_ref = [real(RK)::]

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `gammaIncUpp` with empty input `point` and `kappa`.")

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `gammaIncUpp` with empty input `point` and `kappa`.")

        gammaIncUpp = getGammaIncUpp(point, kappa, tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield an empty `gammaIncUpp` with empty input `point` and `kappa` and `tol`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [EPS]
        point = [0._RK]
        kappa = [1.0_RK]
        GammaIncUpp_ref = [1._RK]

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [0.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [0.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa, tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [0.0_RK], shape = [1._RK], tol = EPS`.")

        block
            use pm_kind, only: LK
            logical(LK) :: Failed(size(point))
            call setGammaIncUppContFracNR(gammaIncUpp, point, kappa, failed, tol)
            call report()
            call test%assert(assertion, SK_"getGammaIncUppContFrac() must yield a correct `gammaIncUpp` with `point = [0.0_RK], shape = [1._RK], tol = EPS`.")
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [EPS]
        point = [1.0_RK]
        kappa = [1.0_RK]
        GammaIncUpp_ref = 1._RK - [0.632120558828557678404476229838539224_RK]

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa, tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [EPS]
        point = [3.0_RK]
        kappa = [5.0_RK]
        GammaIncUpp_ref = 1._RK - [0.184736755476227933713267943730238524_RK]

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa, tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [EPS]
        point = [7.0_RK]
        kappa = [0.5_RK]
        GammaIncUpp_ref = 1._RK - [0.999817189367018164968240003449102634_RK]

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK]`.")

        gammaIncUpp = getGammaIncUpp(point, kappa, tol)
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = [1.0_RK], shape = [1._RK], tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tol = [EPS]
        point = [7.0_RK]
        kappa = [0.5_RK]
        GammaIncUpp_ref = 1._RK - [0.999817189367018164968240003449102634_RK]

        gammaIncUpp = [getGammaIncUpp(point(1), kappa(1))]
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = 1.0_RK, shape = 1._RK`.")

        gammaIncUpp = getGammaIncUpp(point(1), kappa(1))
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = 1.0_RK, shape = 1._RK`.")

        gammaIncUpp = getGammaIncUpp(point(1), kappa(1), tol(1))
        call report()
        call test%assert(assertion, SK_"getGammaIncUpp() must yield a correct `gammaIncUpp` with `point = 1.0_RK, shape = 1._RK, tol = EPS`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(gammaIncUpp - GammaIncUpp_ref)
            assertion = assertion .and. all(abs(diff) <= tol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "point          ", point
                write(test%disp%unit,"(*(g0,:,', '))") "gammaIncUpp    ", gammaIncUpp
                write(test%disp%unit,"(*(g0,:,', '))") "GammaIncUpp_ref", GammaIncUpp_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
