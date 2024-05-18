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
!>  This include file contains the implementation of procedures in [pm_distUnifPar](@ref pm_distUnifPar).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifParLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     Cub_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK < ndim, SK_"@getUnifRecLogPDF(): The condition `0 < ndim` must hold. ndim = "//getStr(ndim))
        logPDF = -ndim * logLenEdge
#elif   Rec_ENABLED
        logPDF = -sum(logLenEdge)
#elif   Par_ENABLED
        integer(IK) :: info
        real(RKC) :: gramian(size(repmat, 1, IK), size(repmat, 2, IK))
        CHECK_ASSERTION(__LINE__, size(repmat, 1, IK) == size(repmat, 2, IK), SK_"@getUnifRecLogPDF(): The condition `size(repmat, 1) == size(repmat, 2)` must hold. shape(repmat) = "//getStr(shape(repmat,IK)))
        gramian = matmul(transpose(repmat), repmat)
        call setMatDetSqrtLog(gramian, uppDia, logPDF, info, gramian, transHerm)
        if (info /= 0_IK) error stop SK_"@getUnifRecLogPDF(): The specified input parallelepiped `repmat` is singular with zero determinant."
        !logPDF = -sum(log(gramian))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifParRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

#if     Cub_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK < ndim, SK_"@getUnifParRand(): The condition `0 < ndim` must hold. ndim = "//getStr(ndim))
#endif
        call setUnifRand(rand)
#if     DU_ENABLED
        call setUnifParRand(rand, ub)
#elif   LU_ENABLED
        call setUnifParRand(rand, lb, ub)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifParRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        ! Define the default lower bound.
#if     DU_ENABLED
        real(RKC), parameter :: lb = 0._RKC
#endif
        ! The input uniform random number must be in range `[0, 1)`.
        CHECK_ASSERTION(__LINE__, all(0._RKC <= rand .and. rand < 1._RKC), SK_"@setUnifParRand(): The condition `all(0. <= rand .and. rand < 1.)` must hold. rand = "//getStr(rand))
        ! Perform checks.
#if     Cub_ENABLED
#define ALL
#elif   Rec_ENABLED || Par_ENABLED
        ! Check the length of `lb` and `ub` against `rand`.
        CHECK_ASSERTION(__LINE__, all(size(rand, 1, IK) == shape(ub, IK)), SK_"@setUnifParRand(): The condition `all(size(rand) == shape(ub))` must hold. size(rand), shape(ub) = "//getStr([size(rand, 1, IK), shape(ub, IK)]))
#if     LU_ENABLED
        CHECK_ASSERTION(__LINE__, size(rand, 1, IK) == size(lb, 1, IK), SK_"@setUnifParRand(): The condition `size(rand) == size(lb)` must hold. size(rand), size(lb) = "//getStr([size(rand, 1, IK), size(lb, 1, IK)]))
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Generate the random vector.
#if     Par_ENABLED
#if     CHECK_ENABLED
        block
            integer(IK) :: idim
            do idim = 1, size(ub, 1, IK)
                CHECK_ASSERTION(__LINE__, 0._RKC < norm2(ub(:,idim)), SK_"@setUnifParRand(): The condition `0. < norm2(ub(:,idim))` must hold. idim, ub = "//getStr(idim)//SK_", "//getStr(ub))
            end do
        end block
#endif
        rand = lb + matmul(ub, rand)
#elif   Rec_ENABLED || Cub_ENABLED
        CHECK_ASSERTION(__LINE__, ALL(lb /= ub), SK_"@setUnifParRand(): The condition `all(lb /= ub)` must hold. lb, ub = "//getStr([lb, ub]))
        rand = (1._RKC - rand) * lb + rand * ub
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  ALL