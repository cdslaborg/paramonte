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
!>  This include file contains implementations of the procedures in module [pm_distNormShell](@ref pm_distNormShell).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHECK_POSITIVE_RADIUS \
CHECK_ASSERTION(__LINE__, all([0._RKG < radius]), \
SK_"@getNormShellLogUDF(): The condition `all([0._RKG < radius])` must hold. width = "//getStr(radius)) ! fpp

#define CHECK_POSITIVE_WIDTH \
CHECK_ASSERTION(__LINE__, all([0._RKG < width]), \
SK_"@getNormShellLogUDF(): The condition `all([0._RKG < width])` must hold. width = "//getStr(width)); ! fpp

#define CHECK_LENGTH_RADIUS \
CHECK_ASSERTION(__LINE__, size(radius, 1, IK) == size(center, 2, IK), \
SK_"@getNormShellLogUDF(): The condition `size(radius, 1) == size(center, 2)` must hold. size(radius, 1), size(center, 2) = "//\
getStr([size(radius, 1, IK), size(center, 2, IK)])) ! fpp

#define CHECK_LENGTH_WIDTH \
CHECK_ASSERTION(__LINE__, size(width, 1, IK) == size(center, 2, IK), \
SK_"@getNormShellLogUDF(): The condition `size(width, 1) == size(center, 2)` must hold. size(width, 1), size(center, 2) = "//\
getStr([size(width, 1, IK), size(center, 2, IK)])) ! fpp

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getNormShellLogUDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: ZERO = 0._RKG, ONE = 1._RKG !, TWO = 2._RKG !, PI = acos(-1._RKG)
#if     CI_ENABLED
        CHECK_ASSERTION(__LINE__, size(X, 1, IK) == size(center, 1, IK), \
        SK_"@getNormShellLogUDF(): The condition `size(X, 1) == size(center, 1)` must hold. size(X, 1), size(center, 1) = "//\
        getStr([size(X, 1, IK), size(center, 1, IK)])) ! fpp

        CHECK_ASSERTION(__LINE__, all(size(X, 1, IK) == [size(invCov, 1, IK), size(invCov, 2, IK)]), \
        SK_"@getNormShellLogUDF(): The condition `all(size(X, 1) == [size(invCov, 1), size(invCov, 2)])` must hold. size(X, 1), shape(invCov) = "//\
        getStr([size(X, 1), shape(invCov)])) ! fpp

        CHECK_ASSERTION(__LINE__, size(invCov, rank(invCov)) == size(center, rank(center)), \
        SK_"@getNormShellLogUDF(): The condition `size(invCov, rank(invCov)) == size(center, rank(invCov))` must hold. shape(invCov), shape(center) = "//\
        getStr([shape(invCov), shape(center)])) ! fpp
#elif   !DD_ENABLED
#error  "Unrecognized interface."
#endif
        ! Compute the density function(s).
#if     D1_ENABLED && (One_ENABLED || (Mix_ENABLED && CI_ENABLED))
#if     DD_ENABLED
#define MAHAL_SQ dot_product(X, X)
#elif   CI_ENABLED
#define MAHAL_SQ getDisMahalSq(X, invCov, center)
#else
#error  "Unrecognized interface."
#endif
        if (present(width)) then
            CHECK_POSITIVE_WIDTH
            if (present(radius)) then
                CHECK_POSITIVE_RADIUS
                logUDF = -.5_RKG * ((sqrt(MAHAL_SQ) - radius) / width)**2
            else
                logUDF = -.5_RKG * ((sqrt(MAHAL_SQ) - ONE) / width)**2
            end if
        else
            if (present(radius)) then
                CHECK_POSITIVE_RADIUS
                logUDF = -.5_RKG *  (sqrt(MAHAL_SQ) - radius)**2
            else
                logUDF = -.5_RKG *  (sqrt(MAHAL_SQ) - ONE)**2
            end if
        end if
#elif   Mix_ENABLED && D1_ENABLED && DD_ENABLED
        CHECK_ASSERTION(__LINE__, size(width, 1, IK) == size(radius, 1, IK), \
        SK_"@getNormShellLogUDF(): The condition `size(width) == size(radius)` must hold. size(width), size(radius) = "//\
        getStr([size(width, 1, IK), size(radius, 1, IK)])) ! fpp
        CHECK_POSITIVE_RADIUS
        CHECK_POSITIVE_WIDTH
        logUDF = -.5_RKG * ((sqrt(sum(X**2)) - radius) / width)**2
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  CHECK_POSITIVE_RADIUS
#undef  CHECK_POSITIVE_WIDTH
#undef  CHECK_LENGTH_RADIUS
#undef  CHECK_LENGTH_WIDTH
#undef  MAHAL_SQ