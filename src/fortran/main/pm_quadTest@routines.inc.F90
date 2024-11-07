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
!>  This include file contains procedure implementation of the generic interfaces of [pm_quadTest](@ref pm_quadTest).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     test_isFailedQuad_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)                 :: numFuncEval
        real(RKG)                   :: integral, abserr, lb, ub, truth
        real(RKG)   , allocatable   :: break(:)

        truth = real(integrand%integral, RKG)
        lb = real(integrand%lb, RKG)
        ub = real(integrand%ub, RKG)

        call disp%show("integrand%desc")
        call disp%show( integrand%desc , deliml = SK_"""" )

        if (present(abstol)) then
            call disp%show("abstol")
            call disp%show( abstol )
        end if

        if (present(reltol)) then
            call disp%show("reltol")
            call disp%show( reltol )
        end if

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("! "//integrand%desc)
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("! Regular adaptive global quadrature.")
        call disp%skip()

        call disp%show("if (isFailedQuad(getFunc, lb, ub, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)")
                        if (isFailedQuad(getFunc, lb, ub, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)
        call disp%show("getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr)")
        call disp%show( getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr) )
        call disp%show("numFuncEval")
        call disp%show( numFuncEval )
        call disp%skip()

        call disp%skip()
        call disp%show("! Assisted adaptive global quadrature by the Wynn Epsilon extrapolation.")
        call disp%skip()

        call disp%show("if (isFailedQuad(getFunc, lb, ub, weps, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)")
                        if (isFailedQuad(getFunc, lb, ub, weps, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)
        call disp%show("getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr)")
        call disp%show( getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr) )
        call disp%show("numFuncEval")
        call disp%show( numFuncEval )
        call disp%skip()

        if (allocated(integrand%break)) then

            call disp%skip()
            call disp%show("! Assisted adaptive global quadrature by specifying points of difficulties of the integrand.")
            call disp%skip()

            call disp%show("break = integrand%break")
                            break = integrand%break
            call disp%show("break")
            call disp%show( break )
            call disp%show("if (isFailedQuad(getFunc, lb, ub, break, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)")
                            if (isFailedQuad(getFunc, lb, ub, break, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)
            call disp%show("getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr)")
            call disp%show( getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr) )
            call disp%show("numFuncEval")
            call disp%show( numFuncEval )
            call disp%skip()

        elseif (allocated(integrand%wcauchy)) then

            call disp%skip()
            call disp%show("! Assisted adaptive global quadrature by specifying Cauchy weight of the integrand.")
            call disp%skip()

            call disp%show("integrand%wcauchy%cs")
            call disp%show( integrand%wcauchy%cs )
            call disp%show("if (isFailedQuad(getFuncUnweighted, lb, ub, integrand%wcauchy, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)")
                            if (isFailedQuad(getFuncUnweighted, lb, ub, integrand%wcauchy, integral, abserr, neval = numFuncEval)) call disp%show('       ******** Integration did NOT converge. ********', tmsize = 1_IK, bmsize = 1_IK)
            call disp%show("getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr)")
            call disp%show( getStr([truth, integral, abserr, abs(integral - truth)])//SK_' (unbiased)? '//getStr(abs(integral - truth) <= abserr) )
            call disp%show("numFuncEval")
            call disp%show( numFuncEval )
            call disp%skip()

        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   test_getQuadErr_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)                 :: nintmax_def
        integer(IK)                 :: err, numFuncEval, numInterval
        real(RKG)                   :: integral, abserr, abstol, reltol, lb, ub, truth
        real(RKG)   , allocatable   :: nodeK(:), weightK(:), weightG(:), sinfo(:,:), break(:)
        integer(IK) , allocatable   :: sindex(:)

        abstol = getOption(0._RKG, atol)
        reltol = getOption(epsilon(0._RKG)**0.66, rtol)
        nintmax_def = getOption(2000_IK, nintmax)
        call setResized(sindex, nintmax_def)
        call setResized(sinfo, [4_IK, nintmax_def])
        truth = real(integrand%integral, RKG)
        lb = real(integrand%lb, RKG)
        ub = real(integrand%ub, RKG)

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("! "//integrand%desc)
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

#define DESC_INTEGRAND \
call disp%show("integrand%desc"); \
call disp%show( integrand%desc , deliml = SK_"""" ); \
call disp%show("[abstol, reltol]"); \
call disp%show( [abstol, reltol] ); \
call disp%show("[lb, ub]"); \
call disp%show( [lb, ub] ); \
call disp%skip();

#define DISP_INTEGRAL \
call disp%show("if (err /= 0_IK) call disp%show(getStrUpper(SK_'        **** integration failed: err = '//getStr(err)//SK_' ****'), tmsize = 1_IK, bmsize = 1_IK)"); \
                if (err /= 0_IK) call disp%show(getStrUpper(SK_'        **** integration failed: err = '//getStr(err)//SK_' ****'), tmsize = 1_IK, bmsize = 1_IK); \
call disp%show("getStr([truth, integral, abserr, abs(integral - truth)])//SK_' >= 1 (unbiased)? '//getStr(abs(integral - truth) <= abserr)"); \
call disp%show( getStr([truth, integral, abserr, abs(integral - truth)])//SK_' >= 1 (unbiased)? '//getStr(abs(integral - truth) <= abserr) ); \
call disp%show("[numFuncEval, numInterval]"); \
call disp%show( [numFuncEval, numInterval] ); \
call disp%skip();

#define DISP_ARB_WEIGHT \
call disp%show("call setNodeWeightGK(order = 35_IK, nodeK = nodeK, weightK = weightK, weightG = weightG) ! Compute the 35-71 Gauss-Kronrod nodes and weights."); \
                call setNodeWeightGK(order = 35_IK, nodeK = nodeK, weightK = weightK, weightG = weightG);

! Define the integration with arbitrary Gauss-Kronrod rule for a given integration range.
#define TRIPLET_ARB(LBV,LBS,UBV,UBS) \
DESC_INTEGRAND \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, nodeK, weightK, weightG, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol, nodeK, weightK, weightG, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
DESC_INTEGRAND \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, nodeK, weightK, weightG, weps, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol, nodeK, weightK, weightG, weps, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
if (allocated(integrand%break)) then; \
    DESC_INTEGRAND \
    call disp%show("break = integrand%break"); \
                    break = integrand%break; \
    call disp%show("break"); \
    call disp%show( break ); \
    DISP_ARB_WEIGHT \
    call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, nodeK, weightK, weightG, break, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                    err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol, nodeK, weightK, weightG, break, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
    DISP_INTEGRAL \
elseif (allocated(integrand%wcauchy)) then; \
    DESC_INTEGRAND \
    call disp%show("integrand%wcauchy%cs"); \
    call disp%show( integrand%wcauchy%cs ); \
    DISP_ARB_WEIGHT \
    call disp%show("err = getQuadErr(getFuncUnweighted, "//LBS//", "//UBS//", abstol, reltol, nodeK, weightK, weightG, integrand%wcauchy, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                    err = getQuadErr(getFuncUnweighted,    LBV,       UBV   , abstol, reltol, nodeK, weightK, weightG, integrand%wcauchy, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
    DISP_INTEGRAL \
end if;

! Define the integration with predefined Gauss-Kronrod rule for a given integration range.
#define TRIPLET_GKX(LBV,LBS,UBV,UBS,GKV,GKS) \
DESC_INTEGRAND \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, "//GKS//", integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol,    GKV   , integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
DESC_INTEGRAND \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, "//GKS//", weps, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol,    GKV   , weps, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
if (allocated(integrand%break)) then; \
DESC_INTEGRAND \
call disp%show("break = integrand%break"); \
                break = integrand%break; \
call disp%show("break"); \
call disp%show( break ); \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFunc, "//LBS//", "//UBS//", abstol, reltol, "//GKS//", break, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFunc,    LBV,       UBV   , abstol, reltol,    GKV   , break, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
elseif (allocated(integrand%wcauchy)) then; \
DESC_INTEGRAND \
call disp%show("integrand%wcauchy%cs"); \
call disp%show( integrand%wcauchy%cs ); \
DISP_ARB_WEIGHT \
call disp%show("err = getQuadErr(getFuncUnweighted, "//LBS//", "//UBS//", abstol, reltol, "//GKS//", integrand%wcauchy, integral, abserr, sinfo, sindex, numFuncEval, numInterval)"); \
                err = getQuadErr(getFuncUnweighted,    LBV,       UBV   , abstol, reltol,    GKV   , integrand%wcauchy, integral, abserr, sinfo, sindex, numFuncEval, numInterval); \
DISP_INTEGRAL \
end if;

        if (getInfNeg(0._RKG) < lb .and. ub < getInfPos(0._RKG)) then

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using arbitrary Gauss-Kronrod nodes and weights (here: 35-71).")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_ARB(lb,"lb",ub,"ub");

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 30-61 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK61,"GK61")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 25-51 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK51,"GK51")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 20-41 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK41,"GK41")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 15-31 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK31,"GK31")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 10-21 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK21,"GK21")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod  7-15 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",ub,"ub",GK15,"GK15")

        elseif (getInfNeg(0._RKG) < lb .and. huge(0._RKG) <= ub) then

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using arbitrary Gauss-Kronrod nodes and weights (here: 35-71).")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_ARB(lb,"lb",pinf,"pinf");

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 30-61 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK61,"GK61")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 25-51 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK51,"GK51")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 20-41 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK41,"GK41")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 15-31 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK31,"GK31")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 10-21 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK21,"GK21")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod  7-15 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(lb,"lb",pinf,"pinf",GK15,"GK15")

        elseif (lb <= -huge(0._RKG) .and. ub < getInfPos(0._RKG)) then

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using arbitrary Gauss-Kronrod nodes and weights (here: 35-71).")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_ARB(ninf,"ninf",ub,"ub");

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 30-61 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK61,"GK61")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 25-51 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK51,"GK51")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 20-41 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK41,"GK41")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 15-31 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK31,"GK31")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 10-21 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK21,"GK21")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod  7-15 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",ub,"ub",GK15,"GK15")

        elseif (lb <= -huge(0._RKG) .and. huge(0._RKG) <= ub) then

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using arbitrary Gauss-Kronrod nodes and weights (here: 35-71).")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_ARB(ninf,"ninf",pinf,"pinf");

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 30-61 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK61,"GK61")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 25-51 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK51,"GK51")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 20-41 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK41,"GK41")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 15-31 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK31,"GK31")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod 10-21 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK21,"GK21")

            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%show("! Using predefined Gauss-Kronrod  7-15 nodes and weights.")
            call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            call disp%skip()

            TRIPLET_GKX(ninf,"ninf",pinf,"pinf",GK15,"GK15")

        end if

#undef  COMPUTE_INTEGRAL
#undef  DISP_ARB_WEIGHT
#undef  DESC_INTEGRAND
#undef  DISP_INTEGRAL
#undef  TRIPLET_GKX
#undef  TRIPLET_ARB

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

    contains

        function getFunc(x) result(func)
            use pm_kind, only: RKH
            real(RKG)    , intent(in)    :: x
            real(RKG)                    :: func
            func = real(integrand%get(real(x, RKH)), RKG)
            if (allocated(integrand%wcauchy)) func = func / (x - real(integrand%wcauchy%cs, RKG))
        end function

        function getFuncUnweighted(x) result(func)
            use pm_kind, only: RKH
            real(RKG)    , intent(in)    :: x
            real(RKG)                    :: func
            func = real(integrand%get(real(x, RKH)), RKG)
        end function
