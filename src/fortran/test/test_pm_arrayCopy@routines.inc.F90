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
!>  This file contains procedure implementations of [test_pm_arrayCopy](@ref test_pm_arrayCopy).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(x) len(x, kind = IK)
#define GET_INDEX(i) i:i
#else
#define GET_SIZE(x) size(x, kind = IK)
#define GET_INDEX(i) i
#endif

#if     SK_ENABLED && D0_ENABLED
#define ALL
        use pm_str, only: getCharSeq, getCharVec
        character(:,SKG), allocatable   :: From, To, To_ref
        character(1,SKG), parameter     :: low = SKG_"a"
        character(1,SKG), parameter     :: upp = SKG_"z"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: From, To, To_ref
        character(2,SKG), parameter                 :: low = SKG_"aa"
        character(2,SKG), parameter                 :: upp = SKG_"zz"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: From, To, To_ref
        integer(IKG)    , parameter                 :: low = -huge(1_IKG)
        integer(IKG)    , parameter                 :: upp = +huge(1_IKG)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: From, To, To_ref
        logical(LKG)    , parameter                 :: low = .false._LKG
        logical(LKG)    , parameter                 :: upp = .true._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: From, To, To_ref
        complex(CKG)    , parameter                 :: low = -cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: upp = +cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: From, To, To_ref
        real(RKG)       , parameter                 :: low = -huge(0._RKG)
        real(RKG)       , parameter                 :: upp = +huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: itest, lenFrom, lenTo, lenCopy

        !%%%%%%%%%%%%%%%%%%%%%
#if     setCopyIndexed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK), allocatable :: indexF(:), indexT(:)

        assertion = .true._LK
        do itest = 1, 100
            lenTo = getUnifRand(1_IK, 10_IK)
            lenFrom = getUnifRand(1_IK, 10_IK)
            lenCopy = getUnifRand(0_IK, max(lenFrom, lenTo))
            indexF = getUnifRand(1_IK, lenFrom, lenCopy)
            indexT = getUnifRand(1_IK, lenTo, lenCopy)
            call setResized(From, lenFrom)
            call setResized(To, lenTo)
            call setUnifRand(From)
            call setUnifRand(To)
            To_ref = To
            call setCopyIndexed_ref()
            call setCopyIndexed(From, To, indexF, indexT)
            call report()
        end do

    contains

        subroutine setCopyIndexed_ref()
            integer(IK) :: i
            do i = 1_IK, lenCopy
                To_ref(GET_INDEX(indexT(i))) = From(GET_INDEX(indexF(i)))
            end do
        end subroutine

        subroutine report()
            use pm_io, only: display_type
            type(display_type) :: disp
            assertion = assertion .and. logical(ALL(To IS_EQUAL To_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("From")
                call disp%show( From )
                call disp%show("indexF")
                call disp%show( indexF )
                call disp%show("indexT")
                call disp%show( indexT )
                call disp%show("To_ref")
                call disp%show( To_ref )
                call disp%show("To")
                call disp%show( To )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The array must be copied correctly.", int(__LINE__, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setCopyStrided_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: incf, inct

        assertion = .true._LK
        do itest = 1, 100
            do
                incf = getUnifRand(-5_IK, 5_IK)
                inct = getUnifRand(-5_IK, 5_IK)
                if (incf == 0_IK .and. inct == 0_IK) cycle
                exit
            end do
            lenCopy = getUnifRand(1_IK, 10_IK)
            lenFrom = merge(lenCopy, (lenCopy - 1_IK) * abs(incf) + 1_IK, incf == 0_IK)
            lenTo   = merge(lenCopy, (lenCopy - 1_IK) * abs(inct) + 1_IK, inct == 0_IK)
            call setResized(From, lenFrom)
            call setResized(To, lenTo)
            call setUnifRand(From)
            call setUnifRand(To)
            To_ref = To

            call setCopyStrided_ref()
            call setCopyStrided(From, To, incf, inct)
            call report()
        end do
        
    contains
        
        subroutine setCopyStrided_ref()
            integer(IK) :: ifrom, ito
            if (incf > 0_IK) then
                ito = 1_IK
                if (inct < 0_IK) ito = GET_SIZE(To_ref)
                do ifrom = 1_IK, GET_SIZE(From), incf
                    To_ref(GET_INDEX(ito)) = From(GET_INDEX(ifrom))
                    ito = ito + inct
                end do
            elseif (incf < 0_IK) then
                ito = 1_IK
                if (inct < 0_IK) ito = GET_SIZE(To_ref)
                do ifrom = GET_SIZE(From), 1_IK, incf
                    To_ref(GET_INDEX(ito)) = From(GET_INDEX(ifrom))
                    ito = ito + inct
                end do
            elseif (inct > 0_IK) then
                do concurrent(ito = 1_IK : GET_SIZE(To_ref) : inct)
                    To_ref(GET_INDEX(ito)) = From(GET_INDEX(1_IK))
                end do
            elseif (inct < 0_IK) then
                do concurrent(ito = GET_SIZE(To_ref) : 1_IK : inct)
                    To_ref(GET_INDEX(ito)) = From(GET_INDEX(1_IK))
                end do
            end if
        end subroutine
        
        subroutine report()
            use pm_io, only: display_type
            type(display_type) :: disp
            assertion = assertion .and. logical(ALL(To IS_EQUAL To_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("From")
                call disp%show( From )
                call disp%show("incf")
                call disp%show( incf )
                call disp%show("inct")
                call disp%show( inct )
                call disp%show("To_ref")
                call disp%show( To_ref )
                call disp%show("To")
                call disp%show( To )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The array must be copied correctly.", int(__LINE__, IK))
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef GET_INDEX
#undef GET_SIZE
#undef IS_EQUAL
#undef ALL