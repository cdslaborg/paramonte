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
!>  This file contains implementations of procedures [pm_distanceEuclid](@ref pm_distanceEuclid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getDisEuclid_ENABLED && XYZ_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: absx, absy, absz, maximum
        absx = abs(x)
        absy = abs(y)
        absz = abs(z)
        maximum = max(absx, absy, absz)
        if(maximum == 0._TKG .or. maximum > huge(maximum)) then
            ! maximum can be zero for max(0,nan,0) adding all three entries together will make sure nan will not disappear.
            distance = absx + absy + absz
        else
            distance = maximum * sqrt((absx / maximum)**2 + (absy / maximum)**2 + (absz / maximum)**2)
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   getDisEuclid_ENABLED
        !%%%%%%%%%%%%%%%%%%%

#if     D1_D1_ENABLED || D1_D2_ENABLED || D2_D1_ENABLED || D2_D2_ENABLED
#define REF ref,
#elif   D1_XX_ENABLED || D2_XX_ENABLED
#define REF
#else
#error  "Unrecognized interface."
#endif
        if (present(method)) then
            select type (method)
            type is (euclid_type)
                call setDisEuclid(distance, point, REF method)
            type is (euclidu_type)
                call setDisEuclid(distance, point, REF method)
            type is (euclidv_type)
                call setDisEuclid(distance, point, REF method)
            type is (euclidsq_type)
                call setDisEuclid(distance, point, REF method)
            class default
                error stop MODULE_NAME//SK_"@getDisEuclid(): Unsupported input value for `method`." ! LCOV_EXCL_LINE
            end select
            return
        end if
        call setDisEuclid(distance, point, REF euclid)
#undef  REF

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisEuclid_ENABLED && (D0_D1_XX_ENABLED || D0_D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, ndim
        ndim = size(point, 1, IK)
#if     D0_D1_XX_ENABLED
#define GET_DIFF(PNT,REF)PNT
#elif   D0_D1_D1_ENABLED
#define GET_DIFF(PNT,REF)(PNT - REF)
        CHECK_ASSERTION(__LINE__, ndim == size(ref, 1, IK), SK_"@setDisEuclid(): The condition `size(point) == size(ref)` must hold. size(point), size(ref) = "//getStr([ndim, size(ref, 1, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     MED_ENABLED
        block
            real(TKG) :: invmax, maximum
            ! First find the maximum.
            maximum = -huge(maximum)
            do idim = 1_IK, ndim
                invmax = abs(GET_DIFF(point(idim),ref(idim))) ! placeholder.
                if (maximum < invmax) maximum = invmax
            end do
            if(maximum == 0._TKG .or. huge(maximum) < maximum) then
                distance = sum(GET_DIFF(point,ref)) ! Ensure propagation of potential nan in the vector.
            else
                distance = 0._TKG
                invmax = 1._TKG / maximum
                do idim = 1_IK, size(point, 1, IK)
                    distance = distance + (GET_DIFF(point(idim),ref(idim)) * invmax)**2
                end do
                distance = maximum * sqrt(distance)
            end if
        end block
#elif   MEU_ENABLED || MEV_ENABLED || MEQ_ENABLED
        distance = 0._TKG
        do idim = 1_IK, ndim
            distance = distance + GET_DIFF(point(idim),ref(idim))**2
        end do
#if     MEU_ENABLED
        distance = sqrt(distance)
#elif   MEV_ENABLED
        block
            real(TKG) :: volUnitBall
            call setVolUnitBall(volUnitBall, ndim)
            distance = volUnitBall * distance**(.5_TKG * ndim)
        end block
#elif   !MEQ_ENABLED
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
#undef  GET_DIFF

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisEuclid_ENABLED && D1_D1_D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iref, ndim
        ndim = size(point, 1, IK)
        CHECK_ASSERTION(__LINE__, size(ref, 2, IK) == size(distance, 1, IK), SK_"@getDisEuclid(): The condition `size(ref, 2) == size(distance)` must hold. shape(ref), size(distance) = "//getStr([shape(ref, IK), size(distance, 1, IK)]))
        do iref = 1_IK, size(ref, 2, IK)
            call setDisEuclid(distance(iref), point, ref(1:ndim, iref), method)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisEuclid_ENABLED && (D1_D2_XX_ENABLED || D1_D2_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipnt, ndim
        ndim = size(point, 1, IK)
        CHECK_ASSERTION(__LINE__, size(point, 2, IK) == size(distance, 1, IK), SK_"@getDisEuclid(): The condition `size(point, 2) == size(distance)` must hold. shape(point), size(distance) = "//getStr([shape(point, IK), size(distance, 1, IK)]))
        do ipnt = 1_IK, size(point, 2, IK)
#if         D1_D2_XX_ENABLED
            call setDisEuclid(distance(ipnt), point(1:ndim, ipnt), method)
#elif       D1_D2_D1_ENABLED
            call setDisEuclid(distance(ipnt), point(1:ndim, ipnt), ref, method)
#else
#error      "Unrecognized interface."
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisEuclid_ENABLED && D2_D1_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iref, nref, npnt
        npnt = size(point, 1, IK)
        nref = size(ref, 1, IK)
        CHECK_ASSERTION(__LINE__, all(shape(distance, IK) == [npnt, nref]), SK_"@getDisEuclid(): The condition `all(shape(distance) == [size(point), size(ref)])` must hold. shape(distance), shape(point), shape(ref) = "//getStr([shape(distance), shape(point), shape(ref)]))
        do iref = 1_IK, nref
#if         MEQ_ENABLED
            distance(1 : npnt, iref) = (point - ref(iref))**2
#elif       MEV_ENABLED
            distance(1 : npnt, iref) = abs(point - ref(iref)) * 2
#elif       MED_ENABLED || MEU_ENABLED
            distance(1 : npnt, iref) = abs(point - ref(iref))
#else
#error      "Unrecognized interface."
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisEuclid_ENABLED && D2_D2_D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iref, npnt, ndim
        npnt = size(point, 2, IK)
        ndim = size(ref, 1, IK)
        CHECK_ASSERTION(__LINE__, size(point, 1, IK) == size(ref, 1, IK), SK_"@getDisEuclid(): The condition `size(point, 1) == size(ref, 1)` must hold. size(point, 1) == size(ref, 1) = "//getStr([size(point, 1, IK) == size(ref, 1, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(distance, IK) == [npnt, size(ref, 2, IK)]), SK_"@getDisEuclid(): The condition `all(shape(distance) == [size(point, 2), size(ref, 2)])` must hold. shape(distance), shape(point), shape(ref) = "//getStr([shape(distance), shape(point), shape(ref)]))
        do iref = 1_IK, size(ref, 2, IK)
            call setDisEuclid(distance(1 : npnt, iref), point, ref(1 : ndim, iref), method)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisMatEuclid_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     FUL_ENABLED
        type(rdpack_type), parameter :: pack = rdpack_type()
        type(uppLowDia_type), parameter :: subset = uppLowDia_type()
#endif
        if (present(method)) then
            select type (method)
            type is (euclidsq_type)
                call setDisMatEuclid(distance, pack, subset, point, method)
            type is (euclidu_type)
                call setDisMatEuclid(distance, pack, subset, point, method)
            type is (euclid_type)
                call setDisMatEuclid(distance, pack, subset, point, method)
            class default
                error stop MODULE_NAME//SK_"@getDisMatEuclid(): Unsupported input value for `method`." ! LCOV_EXCL_LINE
            end select
            return
        end if
        call setDisMatEuclid(distance, pack, subset, point, euclid)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisMatEuclid_ENABLED && RDP_ENABLED && ULD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Construct the full distance matrix.
        integer(IK) :: ndim, npnt
        integer(IK) :: ipnt, jpnt
        ndim = size(point, 1, IK)
        npnt = size(point, 2, IK)
        CHECK_ASSERTION(__LINE__, all(shape(distance, IK) == [npnt, npnt]), \
        SK_"@setDisMatEuclid(): The condition `all(shape(distance) == shape(point))` must hold. shape(distance), shape(point) = "//getStr([shape(distance, IK), shape(point, IK)]))
        do jpnt = 1_IK, npnt
            distance(jpnt, jpnt) = 0._TKG
            do ipnt = jpnt + 1, npnt
                ! construct the lower triangle element and transpose to upper triangle.
                call setDisEuclid(distance(ipnt, jpnt), point(1 : ndim, ipnt), point(1 : ndim, jpnt), method)
                distance(jpnt, ipnt) = distance(ipnt, jpnt)
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setDisMatEuclid_ENABLED && RDP_ENABLED && ULX_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Construct the full distance matrix excluding diagonal.
        integer(IK), parameter :: doff = 1_IK ! diagonal offset.
        integer(IK) :: ndim, npnt
        integer(IK) :: ipnt, jpnt
        ndim = size(point, 1, IK)
        npnt = size(point, 2, IK)
        CHECK_ASSERTION(__LINE__, all(shape(distance, IK) == [npnt - 1_IK, npnt]), \
        SK_"@setDisMatEuclid(): The condition `all(shape(distance) == [size(point, 2) - 1_IK, size(point, 2)])` must hold. shape(distance), shape(point) = "//\
        getStr([shape(distance, IK), shape(point, IK)]))
        do jpnt = 1_IK, npnt
            do ipnt = jpnt + 1, npnt
                ! construct the lower triangle element and transpose to upper triangle.
                call setDisEuclid(distance(ipnt - doff, jpnt), point(1 : ndim, ipnt), point(1 : ndim, jpnt), method)
                distance(jpnt, ipnt) = distance(ipnt - doff, jpnt)
            end do
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#elif   getDisEuclid_ENABLED && (PNT_ENABLED || REF_ENABLED)
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: i
!        real(TKG)   :: invmax
!        real(TKG)   :: maximum
!#if     REF_ENABLED
!        real(TKG)   :: diff(size(point, 1, IK))
!        CHECK_ASSERTION(__LINE__, size(point, 1, IK) == size(ref, 1, IK), SK_"@getDisEuclid(): The condition `size(point) == size(ref)` must hold. size(point), size(ref) = "//getStr([size(point, 1, IK) == size(ref, 1, IK)]))
!#else
!#define DIFF point
!#endif
!        maximum = -huge(maximum)
!        do i = 1_IK, size(point, 1, IK)
!#if         PNT_ENABLED
!            invmax = abs(point(i))
!            if (maximum < invmax) maximum = invmax
!#elif       REF_ENABLED
!            diff(i) = abs(point(i) - ref(i))
!            if (maximum < diff(i)) maximum = diff(i)
!#else
!#error      "Unrecognized interface."
!#endif
!        end do
!        if(maximum == 0._TKG .or. maximum > huge(maximum)) then
!            distance = sum(point) ! Ensure propagation of potential nan in the vector.
!        else
!            distance = 0._TKG
!            invmax = 1._TKG / maximum
!            do i = 1_IK, size(point, 1, IK)
!                distance = distance + (DIFF(i) * invmax)**2
!            end do
!            distance = maximum * sqrt(distance)
!        end if
!
!        !%%%%%%%%%%%%%%%%%%%%%%%
!#elif   getDisMatEuclid_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%
!
!!#if     LFP_ENABLED
!!#define OPTIONAL_ARG
!!#elif   Sym_ENABLED
!!#define OPTIONAL_ARG diag,
!!#elif   REF_ENABLED
!!#define OPTIONAL_ARG ref,
!!#else
!!#error  "Unrecognized interface."
!!#endif
!        if (present(method)) then
!            select type (method)
!            type is (euclidsq_type)
!                call setDisMatEuclid(distance, point, OPTIONAL_ARG euclidsq)
!                return
!            type is (euclidu_type)
!                call setDisMatEuclid(distance, point, OPTIONAL_ARG euclidu)
!                return
!                call setDisMatEuclid(distance, point, OPTIONAL_ARG euclid)
!            end select
!        end if
!        ! type is (euclid_type)
!        call setDisMatEuclid(distance, point, OPTIONAL_ARG euclidu)
!#endif
!
!        !%%%%%%%%%%%%%%%%%%%%%%%
!#elif   setDisMatEuclid_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: idim
!
!        !%%%%%%%%%%%%%%%%%%%%%%%
!#elif   setDisMatEuclid_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%%%
!
!        integer(IK) :: ndim, npnt
!        integer(IK) :: ipnt, jref
!#if     MUF_ENABLED || MEQ_ENABLED
!        integer(IK) :: idim
!#endif
!#if     LFP_ENABLED
!        integer(IK) :: indx
!        indx = 0_IK
!#endif
!        ndim = size(point, 1, IK)
!        npnt = size(point, 2, IK)
!#if     DEF_ENABLED || ULD_ENABLED
!#define REFERENCE point
!#define LBND jref + 1_IK
!#define NREF npnt
!#elif   REF_ENABLED
!#define NREF size(ref, 2, IK)
!#define REFERENCE ref
!#define LBND 1_IK
!        CHECK_ASSERTION(__LINE__, ndim == size(ref, 1, IK), SK_"@setDisMatEuclid(): The condition `size(point, 1) == size(ref, 1)` must hold. size(point, 1), size(ref, 1) = "//getStr([ndim, size(ref, 1, IK)]))
!#else
!#error  "Unrecognized interface."
!#endif
!        ! \warning
!        ! This loops is well-reasoned. Do not change it blindly.
!        do jref = 1_IK, NREF
!#if         Sym_ENABLED
!            distance(jref, jref) = diag
!#endif
!            do ipnt = LBND, npnt
!#if             LFP_ENABLED
!                indx = indx + 1_IK
!#define         INDEX indx
!#else
!#define         INDEX ipnt, jref
!#endif
!#if             MUF_ENABLED || MEQ_ENABLED
!                distance(INDEX) = 0._TKG
!                do idim = 1_IK, ndim
!                    distance(INDEX) = distance(INDEX) + (point(idim, ipnt) - REFERENCE(idim, jref))**2
!                end do
!#if             MUF_ENABLED
!                distance(INDEX) = sqrt(distance(INDEX))
!#endif
!#elif           S_ENABLED
!                distance(INDEX) = getDisEuclid(point(1:ndim, ipnt), REFERENCE(1:ndim, jref))
!#else
!#error          "Unrecognized interface."
!#endif
!#if             Sym_ENABLED
!                distance(jref, ipnt) = distance(ipnt, jref)
!#endif
!            end do
!        end do
!
!#else
!        !%%%%%%%%%%%%%%%%%%%%%%%%
!#error  "Unrecognized interface."
!        !%%%%%%%%%%%%%%%%%%%%%%%%
!#endif
!#undef  OPTIONAL_ARG
!#undef  REFERENCE
!#undef  INDEX
!#undef  NREF
!#undef  LBND
