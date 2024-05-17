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
!>  This include file contains the implementation of procedures in [pm_distanceMahal](@ref pm_distanceMahal).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Check the positive-definiteness of `invCov`.
#define CHECK_POSDEF(INVCOV,ISAM) \
CHECK_ASSERTION(__LINE__, isMatClass(INVCOV, posdefmat), \
SK_"@setDisMahalSq(): The condition `isMatClass(invCov(:,:,isam), posdefmat)` must hold. isam, invCov(:, :, isam) = "//getStr(ISAM)//SK_", "//getStr(INVCOV))
! Check bounds of `invCov`.
#define CHECK_LEN_INVCOV \
CHECK_ASSERTION(__LINE__, all(size(point, 1, IK) == [size(invCov, 1, IK), size(invCov, 2, IK)]), \
SK_"@setDisMahalSq(): The condition `all(size(point, 1) == [size(invCov, 1), size(invCov, 2)])` must hold. size(point), shape(invCov) = "//\
getStr([size(point, 1, IK), shape(invCov, IK)]))
! Check bounds of `invCov`.
#if     InvDef_ENABLED
#define CHECK_LEN_CENTER
#define CHECK_LEN_CENTER_INVCOV
#elif   InvCen_ENABLED
#define CHECK_LEN_CENTER \
CHECK_ASSERTION(__LINE__, size(point, 1) == size(center, 1, IK), \
SK_"@setDisMahalSq(): The condition `size(point, 1) == size(center, 1)` must hold. size(point, 1), size(center) = "//\
getStr([size(point, 1, IK), size(center, 1, IK)]))
! Check bounds of `invCov`.
#define CHECK_LEN_CENTER_INVCOV \
CHECK_ASSERTION(__LINE__, size(center, rank(center), IK) == size(invCov, rank(invCov), IK), \
SK_"@setDisMahalSq(): The condition `size(center, rank(center)) == size(invCov, rank(invCov))` must hold. shape(center), shape(invCov) = "//\
getStr([shape(center, IK), shape(invCov, IK)]))
#else
#error  "Unrecognized interface."
#endif
        !   \bug
        !   ifort 2021.8
        !   error #6401: The attributes of this name conflict with those made accessible by a USE statement.   [SIZE]
        !   It appears ifort cannot digest `size` intrinsic in the `do concurrent` declaration.
        !intrinsic :: size
        ! Set the type and kind.
#if     RK_ENABLED
#define GET_CONJG(point) point
#define TYPE_KIND real(RKC)
        real(RKC)   , parameter :: ZERO = 0._RKC
#elif   CK_ENABLED
#define GET_CONJG(point) conjg(point)
#define TYPE_KIND complex(CKC)
        complex(CKC), parameter :: ZERO = (0._CKC, 0._CKC)
#else
#error  "Unrecognized interface."
#endif
#if     RK_ENABLED && !D0_ENABLED
        integer(IK) :: idim
#elif   !(CK_ENABLED || D0_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%
#if     D0_ENABLED
        !%%%%%%%%%

#if     RK_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKC < invCov, SK_": The condition `0. < invCov` must hold. invCov = "//getStr(invCov))
#elif   !CK_ENABLED
#error  "Unrecognized interface."
#endif
        ! Compute the distance.
#if     InvDef_ENABLED
        mahalSq = GET_CONJG(point) * invCov * point
#elif   InvCen_ENABLED
        mahalSq = GET_CONJG((point - center)) * invCov * (point - center)
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   One_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim
#if     InvDef_ENABLED
#define PNORMED(IDIM) point(IDIM)
#elif   InvCen_ENABLED
        TYPE_KIND :: pnormed(size(point, 1, IK))
        pnormed = point - center
#endif
        ndim = size(point, 1, IK)

        CHECK_POSDEF(invCov, 1)
        CHECK_LEN_CENTER
        CHECK_LEN_INVCOV

        ! Compute the distances.
#if     CK_ENABLED
        mahalSq = dot_product(PNORMED(1:ndim), matmul(invCov, PNORMED(1:ndim))) ! fpp
#elif   RK_ENABLED
        !mahalSq = dot_product(PNORMED , matmul(PNORMED, invCov)) ! fpp
        mahalSq = ZERO
        do idim = 1_IK, ndim
            mahalSq = mahalSq + PNORMED(idim) * dot_product(PNORMED(1:ndim), invCov(1:ndim, idim))
        end do
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   One_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, ndim
#if     InvDef_ENABLED
#define PNORMED(IDIM) point(IDIM, ipnt)
#elif   InvCen_ENABLED
        TYPE_KIND :: pnormed(size(point, 1, IK))
#endif
        ndim = size(point, 1, IK)

        CHECK_POSDEF(invCov, 1)
        CHECK_LEN_CENTER
        CHECK_LEN_INVCOV
#if     setDisMahalSq_ENABLED
        CHECK_ASSERTION(__LINE__, size(point, 2, IK) == size(mahalSq, 1, IK), \
        SK_"@setDisMahalSq(): The condition `size(point, 2) == size(mahalSq, 1)` must hold. size(point, 2), shape(mahalSq, 1) = "//\
        getStr([size(point, 2, IK), size(mahalSq, 1, IK)]))
#endif

        ! Compute the distances.

        do ipnt = 1_IK, size(point, 2, IK)
#if         InvCen_ENABLED
            pnormed(1:ndim) = point(1:ndim, ipnt) - center
#endif
#if         CK_ENABLED
            mahalSq(ipnt) = dot_product(PNORMED(1:ndim), matmul(invCov, PNORMED(1:ndim))) ! fpp
#elif       RK_ENABLED
            !mahalSq(ipnt) = getDisMahalSq(PNORMED(1:ndim), invCov)
            mahalSq(ipnt) = ZERO
            do idim = 1_IK, ndim
                mahalSq(ipnt) = mahalSq(ipnt) + PNORMED(idim) * dot_product(PNORMED(1:ndim), invCov(1:ndim, idim))
            end do
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   Mix_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam, ndim
#if     InvDef_ENABLED
#define PNORMED(IDIM) point(IDIM)
#elif   InvCen_ENABLED
        TYPE_KIND :: pnormed(size(point, 1, IK))
#endif
        ndim = size(point, 1, IK)

        CHECK_LEN_CENTER_INVCOV
        CHECK_LEN_CENTER
        CHECK_LEN_INVCOV
#if     setDisMahalSq_ENABLED
        CHECK_ASSERTION(__LINE__, size(invCov, 3, IK) == size(mahalSq, 1, IK), \
        SK_"@setDisMahalSq(): The condition `size(invCov, 3) == size(mahalSq, 1)` must hold. size(invCov, 3), shape(mahalSq, 1) = "//\
        getStr([size(invCov, 3, IK), size(mahalSq, 1, IK)]))
#endif

        ! Compute the distances.

        do isam = 1_IK, size(invCov, 3, IK)
            CHECK_POSDEF(invCov(:, :, isam), isam)
#if         InvCen_ENABLED
            pnormed(1:ndim) = point(1:ndim) - center(1:ndim, isam)
#endif
#if         CK_ENABLED
            mahalSq(isam) = dot_product(PNORMED(1:ndim), matmul(invCov(:, :, isam), PNORMED(1:ndim))) ! fpp
#elif       RK_ENABLED
            mahalSq(isam) = ZERO
            do idim = 1_IK, ndim
                mahalSq(isam) = mahalSq(isam) + PNORMED(idim) * dot_product(PNORMED(1:ndim), invCov(1:ndim, idim, isam))
            end do
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   Mix_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, isam, nsam, ndim
#if     InvDef_ENABLED
#define PNORMED(IDIM) point(IDIM, ipnt)
#elif   InvCen_ENABLED
        TYPE_KIND :: pnormed(size(point, 1, IK))
#endif
        ndim = size(point, 1, IK)
        nsam = size(invCov, 3, IK)

        CHECK_LEN_CENTER_INVCOV
        CHECK_LEN_CENTER
        CHECK_LEN_INVCOV
#if     setDisMahalSq_ENABLED
        CHECK_ASSERTION(__LINE__, all([size(invCov, 3, IK), size(point, 2, IK)] == shape(mahalSq, IK)), \
        SK_"@setDisMahalSq(): The condition `all([size(invCov, 3), size(point, 2)] == shape(mahalSq))` must hold. size(invCov, 3), size(point, 2), shape(mahalSq) = "//\
        getStr([size(invCov, 3, IK), size(point, 2, IK), shape(mahalSq, IK)]))
#endif

        ! Compute the distances.

        do ipnt = 1_IK, size(point, 2, IK)
            do isam = 1_IK, nsam
#if             InvCen_ENABLED
                pnormed(1:ndim) = point(1:ndim, ipnt) - center(1:ndim, isam)
#endif
#if             CK_ENABLED
                mahalSq(isam, ipnt) = dot_product(PNORMED(1:ndim), matmul(invCov(:, :, isam), PNORMED(1:ndim))) ! fpp
#elif           RK_ENABLED
                mahalSq(isam, ipnt) = ZERO
                do idim = 1_IK, ndim
                    mahalSq(isam, ipnt) = mahalSq(isam, ipnt) + PNORMED(idim) * dot_product(PNORMED(1:ndim), invCov(1:ndim, idim, isam))
                end do
#endif
            end do
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_CENTER_INVCOV
#undef  CHECK_LEN_CENTER
#undef  CHECK_LEN_INVCOV
#undef  CHECK_POSDEF
#undef  TYPE_KIND
#undef  GET_CONJG
#undef  PNORMED
