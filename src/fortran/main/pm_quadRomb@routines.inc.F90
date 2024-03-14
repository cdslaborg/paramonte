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
!>  This include file contains procedure implementation of the generic interfaces of [pm_quadRomb](@ref pm_quadRomb).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     getQuadRomb_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     Clos_ENABLED
        real(RKC)   , parameter     :: FACTOR = 0.5_RKC
        real(RKC)   , parameter     :: SHRINKAGE = 0.25_RKC
        integer(IK) , parameter     :: POINT_INCREMENT_BASE = 2_IK
#elif   Open_ENABLED || PWRL_ENABLED || NEXP_ENABLED || PEXP_ENABLED || LBIS_ENABLED || UBIS_ENABLED
        real(RKC)   , parameter     :: FACTOR = 1._RKC / 3._RKC
        real(RKC)   , parameter     :: SHRINKAGE = 1._RKC / 9._RKC
        integer(IK) , parameter     :: POINT_INCREMENT_BASE = 3_IK
        real(RKC)                   :: resolutionDouble
#if     !Open_ENABLED
        real(RKC)                   :: lbTrans, ubTrans
#endif
#if     LBIS_ENABLED || UBIS_ENABLED
        real(RKC)                   :: inverseOneMinusExp, exp2InverseOneMinusExp ! 1 / (1 - exponent), exponent / (1 - exponent)
#endif
#else
#error  "Unrecognized interface."
#endif
        integer(IK) , parameter     :: NSTEP = ceiling(log(epsilon(0._RKC)) / log(SHRINKAGE)) ! 31_IK: Beyond this integer, the newer QuadRombAbscissa values become identical, which automatically fail `setInterp()`.
        integer(IK)                 :: refinementStage, km
        real(RKC)                   :: QuadRombAbscissa(NSTEP + 1), QuadRombProxy(NSTEP + 1)
        real(RKC)                   :: integrationRange, nevalNewInverse, resolution, sumFunc, x
        integer(IK)                 :: nevalNew, ieval
#if     EM_ENABLED
        real(RKC)                   :: relerr
#endif
#if     NP_ENABLED && Clos_ENABLED
        neval = 2_IK
#elif   NP_ENABLED
        neval = 1_IK
#endif
        ! Set the final normalization of the integral, only necessary when there is singularity at a limit.
#if     LBIS_ENABLED || UBIS_ENABLED
#define NORMALIZE_QUAD_ROMB quadRomb = quadRomb * inverseOneMinusExp;
#else
#define NORMALIZE_QUAD_ROMB
#endif
        ! Set the transformed limits
#if     Clos_ENABLED || Open_ENABLED
#define GET_FUNC(x)getFunc(x)
#define lbTrans lb
#define ubTrans ub
#elif   PWRL_ENABLED
#define GET_FUNC(x)getFunc(1._RKC / x) / x**2
        lbTrans = 1._RKC / ub
        ubTrans = 1._RKC / lb
#elif   NEXP_ENABLED
#define GET_FUNC(x)getFunc(-log(x)) / x
        lbTrans = exp(-ub)
        ubTrans = exp(-lb)
#elif   PEXP_ENABLED
#define GET_FUNC(x)getFunc(+log(x)) / x
        lbTrans = exp(lb)
        ubTrans = exp(ub)
#elif   LBIS_ENABLED || UBIS_ENABLED
        lbTrans = 0._RKC
        ubTrans = (ub - lb)**(1._RKC + real(Interval%exponent, RKC))
        inverseOneMinusExp = 1._RKC / (1._RKC + real(Interval%exponent, RKC))
        exp2InverseOneMinusExp = -real(Interval%exponent, RKC) * inverseOneMinusExp
        CHECK_ASSERTION(__LINE__, -1._RKC < Interval%exponent .and. Interval%exponent <= 0._RKC, \
        SK_"The conditions `-1. < Interval%exponent .and. Interval%exponent <= 0.` must hold: Interval%exponent = "//getStr(Interval%exponent))
#if     LBIS_ENABLED
#define GET_FUNC(x)x**exp2InverseOneMinusExp * getFunc(lb + x**inverseOneMinusExp)
#elif   UBIS_ENABLED
#define GET_FUNC(x)x**exp2InverseOneMinusExp * getFunc(ub - x**inverseOneMinusExp)
#endif
#endif
        integrationRange = ubTrans - lbTrans
        CHECK_ASSERTION(__LINE__, integrationRange >= 0._RKC, SK_"The input lower and upper integration bounds [lb, ub] must satisfy `lb < ub`.")
        CHECK_ASSERTION(__LINE__, nref > 0_IK, SK_"The input refinement count in the Romberg method must satisfy `nref > 0`. nref = "//getStr(nref))
        CHECK_ASSERTION(__LINE__, nref <= NSTEP, SK_"The input refinement count in the Romberg method must satisfy `nref > 0`. nref, NSTEP = "//getStr([nref, NSTEP]))
        km = nref - 1_IK
        refinementStage = 1_IK
        QuadRombAbscissa(refinementStage) = 1._RKC
#if     Clos_ENABLED
        QuadRombProxy(refinementStage) = 0.5_RKC * integrationRange * (GET_FUNC(lbTrans) + GET_FUNC(ubTrans))
#else
        QuadRombProxy(refinementStage) = integrationRange * GET_FUNC((0.5_RKC * (lbTrans + ubTrans))) ! \warning extra parentheses are important
#endif
        ! Define the evaluation instruction.
#define EVAL_QUAD_ROMB \
if (refinementStage >= nref) then; \
call setExtrap(monopol, QuadRombAbscissa(refinementStage - km : refinementStage), QuadRombProxy(refinementStage - km : refinementStage), 0._RKC, quadRomb, relerr); \
relerr = abs(relerr); \
if (relerr <= tol * abs(quadRomb)) then; \
NORMALIZE_QUAD_ROMB \
return; \
end if; \
end if; \
QuadRombProxy(refinementStage + 1_IK) = QuadRombProxy(refinementStage); \
QuadRombAbscissa(refinementStage + 1_IK) = SHRINKAGE * QuadRombAbscissa(refinementStage);
        ! Compute the integral.
        EVAL_QUAD_ROMB
        do refinementStage = 2_IK, NSTEP
            nevalNew = POINT_INCREMENT_BASE**(refinementStage - 2_IK)
            nevalNewInverse = 1._RKC / real(nevalNew, RKC)
#if         Clos_ENABLED
            resolution = integrationRange * nevalNewInverse
#else
            resolution = integrationRange * nevalNewInverse * FACTOR
            resolutionDouble = resolution + resolution
#endif
            x = lbTrans + 0.5_RKC * resolution
            sumFunc = 0._RKC
            do ieval = 1_IK, nevalNew
                sumFunc = sumFunc + GET_FUNC(x)
#if             !Clos_ENABLED
                x = x + resolutionDouble
                sumFunc = sumFunc + GET_FUNC(x)
#endif
                x = x + resolution
            end do
            QuadRombProxy(refinementStage) = FACTOR * (QuadRombProxy(refinementStage) + integrationRange * sumFunc * nevalNewInverse)
#if         NP_ENABLED && Clos_ENABLED
            neval = neval + nevalNew
#elif       NP_ENABLED
            neval = neval + nevalNew * 2
#endif
            EVAL_QUAD_ROMB
        end do
#if     EM_ENABLED
        error stop __FILE__//"@getQuadRomb(): The Romberg integration failed to converge."
#elif   EP_ENABLED
        relerr = -huge(relerr)
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  NORMALIZE_QUAD_ROMB
#undef  GET_FUNC
#undef  lbTrans
#undef  ubTrans
