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
!>  This file contains the implementations of the tests of module [pm_mathRoot](@ref pm_mathRoot).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKC)
#elif   RK_ENABLED
#define TYPE_KIND real(TKC)
#else
#error  "Unrecognized interface."
#endif
        ! Select root-finding method.
#if     Brent_ENABLED || Def_ENABLED
#define METHOD brent_type
#elif   False_ENABLED
#define METHOD false_type
#elif   Secant_ENABLED
#define METHOD secant_type
#elif   Newton_ENABLED
#define METHOD newton_type
#elif   Halley_ENABLED
#define METHOD halley_type
#elif   Ridders_ENABLED
#define METHOD ridders_type
#elif   TOMS748_ENABLED
#define METHOD toms748_type
#elif   Schroder_ENABLED
#define METHOD schroder_type
#elif   Bisection_ENABLED
#define METHOD bisection_type
#else
#error  "Unrecognized interface."
#endif
        ! Set the default `getFunc()`.
#if     Newton_ENABLED || Halley_ENABLED || Schroder_ENABLED
#define GETFUNC getFuncDiff
#endif
        type(METHOD), parameter :: method = METHOD()
        type(display_type) :: disp
        integer(IK) :: itry, neval, niter
        real(TKC), parameter :: abstol_ref = epsilon(0._TKC)**.8
        TYPE_KIND :: roots(2), lb, ub, xmid
        assertion = .true._LK

        do itry = 1, 100

            roots = getChoice(getRemoved(getRange(-5, 5), 0), 2_IK, unique = .true._LK)
            lb = getChoice([minval(roots, 1), minval(roots, 1) - 1])
            ub = getChoice([maxval(roots, 1), maxval(roots, 1) + 1])
            if (roots .allinrange. [lb, ub]) then
                xmid = .5_TKC * sum(roots(1:2))
                if (getUnifRand()) then
                    lb = xmid
                else
                    ub = xmid
                end if
            end if
            niter = getUnifRand(100, 300)

#if         getRoot_ENABLED
            call testWith(method)
            call testWith(method, abstol_ref)
            call testWith(method, abstol_ref, neval)
            call testWith(method, abstol_ref, neval, niter)
            call testWith(method, abstol_ref, niter = niter)
            call testWith(method, neval = neval, niter = niter)
            call testWith(method, neval = neval)
            call testWith(method, niter = niter)
            call testWith()
            call testWith(abstol = abstol_ref)
            call testWith(abstol = abstol_ref, neval = neval)
            call testWith(abstol = abstol_ref, neval = neval, niter = niter)
            call testWith(abstol = abstol_ref, niter = niter)
            call testWith(neval = neval, niter = niter)
            call testWith(neval = neval)
            call testWith(niter = niter)
#elif       setRoot_ENABLED
            block
                TYPE_KIND :: root, lf, uf
                lf = getFunc(lb)
                uf = getFunc(ub)
                root = getUnifRand(lb - 1, ub + 1) ! this initialization is required in testing of the Newton method.
                call setRoot(method, GETFUNC, root, lb, ub, lf, uf, abstol_ref, neval)
                assertion = assertion .and. (any(abs(roots - root) < abstol_ref * 1000) .or. getOption(0_IK, neval) < 0_IK)
                call report(__LINE__, root, method, abstol_ref * 1000, neval)

                call setRoot(method, GETFUNC, root, lb, ub, lf, uf, abstol_ref, neval, niter)
                assertion = assertion .and. (any(abs(roots - root) < abstol_ref * 1000) .or. getOption(0_IK, neval) < 0_IK)
                call report(__LINE__, root, method, abstol_ref * 1000, neval, niter)
            end block
#else
#error      "Unrecognized interface."
#endif

        end do

    contains

        pure function getFunc(x) result(func)
            real(TKC), intent(in) :: x
            real(TKC) :: func
            func = x**2 - sum(roots(1:2)) * x + product(roots) ! (x - roots(1)) * (x - roots(2))
        end function

#if     Newton_ENABLED || Halley_ENABLED || Schroder_ENABLED
        pure function getFuncDiff(x, order) result(func)
            integer(IK), intent(in) :: order
            real(TKC), intent(in) :: x
            real(TKC) :: func
            if (order == 0) func = getFunc(x)
            if (order == 1) func = 2._TKC * x - sum(roots(1:2))
            if (order == 2) func = 2._TKC
        end function
#endif

#if     getRoot_ENABLED
        subroutine testWith(method, abstol, neval, niter)
            type(METHOD), intent(in), optional :: method
            integer(IK), intent(inout), optional :: neval, niter
            real(TKC), intent(in), optional :: abstol
            real(TKC) :: abstol_def
            integer(IK) :: neval_def
            TYPE_KIND :: root
            if (present(abstol)) then
                abstol_def = abstol * 1000._TKC
            else
                abstol_def = abstol_ref * (abs(lb) + abs(ub)) * 1000._TKC
            end if
            if (present(method)) then
                root = getRoot(method, GETFUNC, lb, ub, abstol, neval, niter)
                assertion = assertion .and. (any(abs(roots - root) < abstol_def) .or. getOption(0_IK, neval) < 0_IK)
                call report(__LINE__, root, method, abstol, neval, niter)
#if             Newton_ENABLED || Halley_ENABLED || Schroder_ENABLED
                !call disp%show("[real(TKC) :: lb, ub, getFunc(lb), getFunc(ub), abstol_def]")
                !call disp%show( [real(TKC) :: lb, ub, getFunc(lb), getFunc(ub), abstol_def] )
                !call disp%show("[roots, getFunc(roots(1)), getFunc(roots(2))]")
                !call disp%show( [roots, getFunc(roots(1)), getFunc(roots(2))] )
                root = getRoot(method, GETFUNC, lb, ub, abstol, neval, niter, init = getUnifRand(lb - 1, ub + 1))
                assertion = assertion .and. (any(abs(roots - root) < abstol_def) .or. getOption(0_IK, neval) < 0_IK)
                call report(__LINE__, root, method, abstol, neval, niter)
#endif
            else
                !print *, "size(roots)", size(roots)
                !print *, "roots", roots
                !print *, "froots", getFunc(roots(0)), getFunc(roots(1)), getFunc(roots(2))
                !print *, "lb, ub", lb, ub
                root = getRoot(getFunc, lb, ub, abstol, neval, niter)
                assertion = assertion .and. (any(abs(roots - root) < abstol_def) .or. getOption(0_IK, neval) < 0_IK)
                call report(__LINE__, root, method, abstol, neval, niter)
            end if
            root = getRoot(getFunc, lb, ub, abstol, neval_def, 1_IK)
            assertion = assertion .and. (any(abs(roots - root) < abstol_def) .or. neval_def <= 0_IK)
            call report(__LINE__, root, method, abstol, neval_def, niter)
        end subroutine
#endif

        subroutine report(line, root, method, abstol, neval, niter)
            type(METHOD), intent(in), optional :: method
            integer(IK), intent(in), optional :: neval, niter
            real(TKC), intent(in), optional :: abstol
            TYPE_KIND, intent(in) :: root
            integer, intent(in) :: line
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("roots")
                call disp%show( roots )
                call disp%show("[lb, root, ub]")
                call disp%show( [lb, root, ub] )
                call disp%show("getFunc(root)")
                call disp%show( getFunc(root) )
                call disp%show("abs(root - roots)")
                call disp%show( abs(root - roots) )
                call disp%show("present(method)")
                call disp%show( present(method) )
                call disp%show("present(abstol)")
                call disp%show( present(abstol) )
                if (present(abstol)) then
                    call disp%show("abstol")
                    call disp%show( abstol )
                end if
                call disp%show("present(neval)")
                call disp%show( present(neval) )
                if (present(neval)) then
                    call disp%show("neval")
                    call disp%show( neval )
                end if
                call disp%show("present(niter)")
                call disp%show( present(niter) )
                if (present(niter)) then
                    call disp%show("niter")
                    call disp%show( niter )
                end if
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"getRoot() must compute a root of the target function correctly.", int(line, IK))
        end subroutine

#undef  GETFUNC
#undef  METHOD