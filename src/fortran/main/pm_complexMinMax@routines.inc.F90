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
!>  This file contains the implementation details of the routines under the
!>  generic interfaces of module [pm_complexAbs](@ref pm_complexAbs).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the initial values and comparison operations.
#if     minval_ENABLED || minloc_ENABLED
#define INIT+huge(0._CKC)
#define COMPARES_WITH>
#define GETLOC minloc
#define GETVAL minval
#elif   maxval_ENABLED || maxloc_ENABLED
#define INIT-huge(0._CKC)
#define COMPARES_WITH<
#define GETLOC maxloc
#define GETVAL maxval
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     min_ENABLED && D0_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        val%re = min(a1%re, a2%re)
        val%im = min(a1%im, a2%im)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   max_ENABLED && D0_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        val%re = max(a1%re, a2%re)
        val%im = max(a1%im, a2%im)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minval_ENABLED || maxval_ENABLED) && D1_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
        val%re = INIT
        val%im = INIT
        do iell = 1, size(array, 1, IK)
            if (val%re COMPARES_WITH array(iell)%re) val%re = array(iell)%re
            if (val%im COMPARES_WITH array(iell)%im) val%im = array(iell)%im
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minval_ENABLED || maxval_ENABLED) && D2_ENABLED && ALL_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        complex(CKC), pointer :: parray(:)
        parray(1 : size(array, kind = IK)) => array
        val = GETVAL(parray)
        nullify(parray)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minval_ENABLED || maxval_ENABLED) && D2_ENABLED && DIM_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim, nsam
        CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(array), SK_": The condition `1 <= dim .and. dim <= rank(array)` must hold. dim, rank(array) = "//getStr([dim, int(rank(array), IK)]))
        ndim = size(array, 3 - dim, IK)
        nsam = size(array, dim, IK)
        if (dim == 1_IK) then
            do concurrent(idim = 1 : ndim)
                val(idim) = GETVAL(array(1 : nsam, idim))
            end do
        elseif (dim == 2_IK) then
            do concurrent(idim = 1 : ndim)
                val(idim) = GETVAL(array(idim, 1 : nsam))
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minloc_ENABLED || maxloc_ENABLED) && D1_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        complex(CKC) :: val
        integer(IK) :: iell
        loc = 0_IK
        val%re = INIT
        val%im = INIT
        do iell = 1, size(array, 1, IK)
            if (val%re COMPARES_WITH array(iell)%re) then
                val%re = array(iell)%re
                loc(1) = iell
            end if
            if (val%im COMPARES_WITH array(iell)%im) then
                val%im = array(iell)%im
                loc(2) = iell
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minloc_ENABLED || maxloc_ENABLED) && D2_ENABLED && ALL_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        complex(CKC), pointer :: parray(:)
        parray(1 : size(array, kind = IK)) => array
        loc = GETLOC(parray)
        nullify(parray)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (minloc_ENABLED || maxloc_ENABLED) && D2_ENABLED && DIM_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim, nsam
        CHECK_ASSERTION(__LINE__, 1 <= dim .and. dim <= rank(array), SK_": The condition `1 <= dim .and. dim <= rank(array)` must hold. dim, rank(array) = "//getStr([dim, int(rank(array), IK)]))
        ndim = size(array, 3 - dim, IK)
        nsam = size(array, dim, IK)
        if (dim == 1_IK) then
            do concurrent(idim = 1 : ndim)
                loc(1 : 2, idim) = GETLOC(array(1 : nsam, idim))
            end do
        elseif (dim == 2_IK) then
            do concurrent(idim = 1 : ndim)
                loc(1 : 2, idim) = GETLOC(array(idim, 1 : nsam))
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  COMPARES_WITH
#undef  GETLOC
#undef  GETVAL
#undef  INIT