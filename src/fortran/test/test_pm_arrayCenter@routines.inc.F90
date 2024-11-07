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
!>  This module contains implementations of the tests of the procedures under the generic interfaces [pm_arrayCenter](@ref pm_arrayCenter).
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
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) len(Array, kind = IK)
#define GET_SIZE(Array) len(Array, kind = IK)
#define ALLOC_ARRAY(Array,lb,ub) allocate(character(ub-lb+1,SKG) :: Array)
#define GEN_UBOLD(ub) ub - lbold + 1_IK
#define GEN_UBNEW(ub) ub - lbnew + 1_IK + getOption(0_IK,lmsize) + getOption(0_IK,rmsize)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#elif   D1_ENABLED
#define ALLOC_ARRAY(Array,lb,ub) allocate(Array(lb:ub))
#define GET_LBOUND(Array) lbound(Array, dim = 1, kind = IK)
#define GEN_UBOUND(Array) ubound(Array, dim = 1, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_UBOLD(ub) ub
#define GEN_UBNEW(ub) ub + getOption(0_IK,lmsize) + getOption(0_IK,rmsize)
#define GEN_LBOLD(lb) lb
#define GEN_LBNEW(lb) lb
#else
#error  "Unrecognized interface."
#endif

#if     SK_ENABLED && D0_ENABLED
#define ALL
        character(:,SKG), allocatable   :: Array, ArrayCentered, ArrayCentered_ref
        character(1,SKG), parameter     :: fill = SKG_"*"
        character(1,SKG), parameter     :: lmfill = SKG_"-"
        character(1,SKG), parameter     :: rmfill = SKG_"+"
        character(1,SKG), parameter     :: fill_def = fill
        character(1,SKG), parameter     :: lmfill_def = lmfill
        character(1,SKG), parameter     :: rmfill_def = rmfill
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, ArrayCentered, ArrayCentered_ref
        character(2,SKG), parameter                 :: fill = SKG_"**"
        character(2,SKG), parameter                 :: lmfill = SKG_"--"
        character(2,SKG), parameter                 :: rmfill = SKG_"++"
        character(2,SKG), parameter                 :: fill_def = fill
        character(2,SKG), parameter                 :: lmfill_def = lmfill
        character(2,SKG), parameter                 :: rmfill_def = rmfill
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, ArrayCentered, ArrayCentered_ref
        logical(LKG)    , parameter                 :: fill = .false._LKG
        logical(LKG)    , parameter                 :: lmfill = .false._LKG
        logical(LKG)    , parameter                 :: rmfill = .false._LKG
        logical(LKG)    , parameter                 :: fill_def = fill
        logical(LKG)    , parameter                 :: lmfill_def = lmfill
        logical(LKG)    , parameter                 :: rmfill_def = rmfill
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, ArrayCentered, ArrayCentered_ref
        integer(IKG)    , parameter                 :: fill = huge(1_IKG)
        integer(IKG)    , parameter                 :: lmfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: rmfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: fill_def = fill
        integer(IKG)    , parameter                 :: lmfill_def = lmfill
        integer(IKG)    , parameter                 :: rmfill_def = rmfill
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, ArrayCentered, ArrayCentered_ref
        complex(CKG)    , parameter                 :: fill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: lmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: rmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: fill_def = fill
        complex(CKG)    , parameter                 :: lmfill_def = lmfill
        complex(CKG)    , parameter                 :: rmfill_def = rmfill
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, ArrayCentered, ArrayCentered_ref
        real(RKG)       , parameter                 :: fill = huge(0._RKG)
        real(RKG)       , parameter                 :: lmfill = huge(0._RKG)
        real(RKG)       , parameter                 :: rmfill = huge(0._RKG)
        real(RKG)       , parameter                 :: fill_def = fill
        real(RKG)       , parameter                 :: lmfill_def = lmfill
        real(RKG)       , parameter                 :: rmfill_def = rmfill
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: cbeg
        integer(IK) :: cend
        integer(IK) :: sizeold
        integer(IK) :: sizenew
        integer(IK) :: sizecentered
        integer(IK) :: lbold, ubold
        integer(IK) :: lbnew, ubnew
        logical(LK) :: menabled

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(fill = fill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, fill = fill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, lmfill = lmfill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, lmfill = lmfill, fill = fill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, rmfill = rmfill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, rmfill = rmfill, fill = fill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, lmfill = lmfill, rmfill = rmfill)
        call runTestsWith(lmsize = 2_IK, rmsize = 3_IK, lmfill = lmfill, rmfill = rmfill, fill = fill)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith (lmsize, rmsize, fill, lmfill, rmfill)

            integer(IK)     , intent(in), optional  :: lmsize, rmsize
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in), optional  :: fill, lmfill, rmfill
#elif       SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in), optional  :: fill, lmfill, rmfill
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in), optional  :: fill, lmfill, rmfill
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in), optional  :: fill, lmfill, rmfill
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in), optional  :: fill, lmfill, rmfill
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in), optional  :: fill, lmfill, rmfill
#else
#error      "Unrecognized interface."
#endif

            if (present(lmsize) .neqv. present(rmsize)) error stop "Internal ParaMonte Testing error occurred: present(lmsize) .neqv. present(rmsize) must be .true."
            menabled = present(lmsize) .and. present(rmsize)

            !>  \bug
            !>  GNU Fortran 10.3 cannot concatenate empty character array of length 2 with a non-empty character array of the same length.
            !>  Fortran runtime error: Different CHARACTER lengths (0/2) in array constructor
            if (present(lmsize) .and. present(rmsize)) then
                if (lmsize == 0_IK .and. rmsize == 0_IK) error stop "Internal ParaMonte Testing error occurred: GNU bug exception."
            end if

            assertion = .true._LK

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Enlarge and center array
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

            lbold = GEN_LBOLD(1_IK)
            ubold = GEN_UBOLD(3_IK)
            lbnew = GEN_LBNEW(1_IK)
            ubnew = GEN_UBNEW(4_IK)
            sizeold = ubold - lbold + 1_IK
            sizenew = ubnew - lbnew + 1_IK
            sizecentered = sizenew - getOption(0_IK, lmsize) - getOption(0_IK, rmsize)
            ALLOC_ARRAY(Array,lbold,ubold)
            ALLOC_ARRAY(ArrayCentered,lbnew,ubnew)
            ALLOC_ARRAY(ArrayCentered_ref,lbnew,ubnew)
#if         SK_ENABLED && D0_ENABLED
            call setUnifRand(Array, repeat(SKG_"A", len(Array)), repeat(SKG_"Z",len(Array)))
            ArrayCentered_ref = getRepeat(lmfill, lmsize)//Array//fill_def//getRepeat(rmfill, rmsize)
#elif       SK_ENABLED && D1_ENABLED
            Array(:) = [SKG_"AA", SKG_"BB", SKG_"CC"]
            if (getOption(0_IK,lmsize) > 0_IK) then
                ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), Array, fill_def, getRepeat(rmfill,rmsize)]
            else
                ArrayCentered_ref(:) = [Array, fill_def]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array(:) = [.true._LKG, .true._LKG, .true._LKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       IK_ENABLED && D1_ENABLED
            Array(:) = [1_IKG, 2_IKG, 3_IKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       CK_ENABLED && D1_ENABLED
            Array(:) = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG)]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       RK_ENABLED && D1_ENABLED
            Array(:) = [1._RKG, 2._RKG, 3._RKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), Array, fill_def, getRepeat(rmfill,rmsize)]
#endif
            cbeg = lbold + getOption(0_IK,lmsize) + (sizecentered - sizeold) / 2_IK
            cend = cbeg + sizeold - 1_IK

#if         setCentered_ENABLED
            if (menabled) then
                call setCentered(ArrayCentered, Array, lmsize, rmsize, fill, lmfill, rmfill)
            else
                call setCentered(ArrayCentered, Array, fill)
            end if

            assertion = assertion .and. GET_LBOUND(ArrayCentered) == lbnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the lower bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. GEN_UBOUND(ArrayCentered) == ubnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))
#elif       getCentered_ENABLED
            if (menabled) then
                ArrayCentered = getCentered(Array, sizecentered, lmsize, rmsize, fill, lmfill, rmfill)
            else
                ArrayCentered = getCentered(Array, sizecentered, fill)
            end if
#else
#error      "Unrecognized interface."
#endif

            assertion = assertion .and. GET_SIZE(ArrayCentered) == sizenew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must yield an array of proper size, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

#if         setCentered_ENABLED
            assertion = assertion .and. ALL(ArrayCentered(cbeg:cend) IS_EQUAL ArrayCentered_ref(cbeg:cend))
#elif       getCentered_ENABLED
            assertion = assertion .and. ALL(ArrayCentered(cbeg-lbold+1_IK:cend-lbold+1_IK) IS_EQUAL ArrayCentered_ref(cbeg:cend))
#endif
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            if (present(fill)) then
            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly fill the new elements with `fill` = "//getStr(present(fill)), int(__LINE__, IK))
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Enlarge and center array perfectly
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

            lbold = GEN_LBOLD(1_IK)
            ubold = GEN_UBOLD(3_IK)
            lbnew = GEN_LBNEW(1_IK)
            ubnew = GEN_UBNEW(5_IK)
            sizeold = ubold - lbold + 1_IK
            sizenew = ubnew - lbnew + 1_IK
            sizecentered = sizenew - getOption(0_IK,lmsize) - getOption(0_IK,rmsize)
            ALLOC_ARRAY(Array,lbold,ubold)
            ALLOC_ARRAY(ArrayCentered,lbnew,ubnew)
            ALLOC_ARRAY(ArrayCentered_ref,lbnew,ubnew)
#if         SK_ENABLED && D0_ENABLED
            call setUnifRand(Array, repeat(SKG_"A",len(Array)), repeat(SKG_"Z",len(Array)))
            ArrayCentered_ref = getRepeat(lmfill,lmsize)//fill_def//Array//fill_def//getRepeat(rmfill,rmsize)
#elif       SK_ENABLED && D1_ENABLED
            Array(:) = [SKG_"AA", SKG_"BB", SKG_"CC"]
            if (getOption(0_IK,lmsize) > 0_IK) then
                ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), fill_def, Array, fill_def, getRepeat(rmfill,rmsize)]
            else
                ArrayCentered_ref(:) = [fill_def, Array, fill_def]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array(:) = [.true._LKG, .true._LKG, .true._LKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), fill_def, Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       IK_ENABLED && D1_ENABLED
            Array(:) = [1_IKG, 2_IKG, 3_IKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), fill_def, Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       CK_ENABLED && D1_ENABLED
            Array(:) = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG)]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), fill_def, Array, fill_def, getRepeat(rmfill,rmsize)]
#elif       RK_ENABLED && D1_ENABLED
            Array(:) = [1._RKG, 2._RKG, 3._RKG]
            ArrayCentered_ref(:) = [getRepeat(lmfill,lmsize), fill_def, Array, fill_def, getRepeat(rmfill,rmsize)]
#endif
            cbeg = lbold + getOption(0_IK,lmsize) + (sizecentered - sizeold) / 2_IK
            cend = cbeg + sizeold - 1_IK

#if         setCentered_ENABLED
            if (menabled) then
                call setCentered(ArrayCentered, Array, lmsize, rmsize, fill, lmfill, rmfill)
            else
                call setCentered(ArrayCentered, Array, fill)
            end if

            assertion = assertion .and. GET_LBOUND(ArrayCentered) == lbnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the lower bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. GEN_UBOUND(ArrayCentered) == ubnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))
#elif       getCentered_ENABLED
            if (menabled) then
                ArrayCentered = getCentered(Array, sizecentered, lmsize, rmsize, fill, lmfill, rmfill)
            else
                ArrayCentered = getCentered(Array, sizecentered, fill)
            end if
#else
#error      "Unrecognized interface."
#endif

            assertion = assertion .and. GET_SIZE(ArrayCentered) == sizenew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must yield an array of proper size, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

#if         setCentered_ENABLED
            assertion = assertion .and. ALL(ArrayCentered(cbeg:cend) IS_EQUAL ArrayCentered_ref(cbeg:cend))
#elif       getCentered_ENABLED
            assertion = assertion .and. ALL(ArrayCentered(cbeg-lbold+1_IK:cend-lbold+1_IK) IS_EQUAL ArrayCentered_ref(cbeg:cend))
#endif
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            if (present(fill)) then
            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly fill the new elements with `fill` = "//getStr(present(fill)), int(__LINE__, IK))
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Shrink and center array
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

            lbold = GEN_LBOLD(1_IK)
            ubold = GEN_UBOLD(4_IK)
            lbnew = GEN_LBNEW(1_IK)
            ubnew = GEN_UBNEW(3_IK)
            sizeold = ubold - lbold + 1_IK
            sizenew = ubnew - lbnew + 1_IK
            sizecentered = sizenew - getOption(0_IK,lmsize) - getOption(0_IK,rmsize)
            ALLOC_ARRAY(Array,lbold,ubold)
            ALLOC_ARRAY(ArrayCentered,lbnew,ubnew)
            ALLOC_ARRAY(ArrayCentered_ref,lbnew,ubnew)
#if         SK_ENABLED && D0_ENABLED
            call setUnifRand(Array, repeat(SKG_"A",len(Array)), repeat(SKG_"Z",len(Array)))
            ArrayCentered_ref(lbnew:ubnew) = getRepeat(lmfill,lmsize)//Array(lbold:ubold-1)//getRepeat(rmfill,rmsize)
#elif       SK_ENABLED && D1_ENABLED
            Array(:) = [SKG_"AA", SKG_"BB", SKG_"CC", SKG_"DD"]
            if (getOption(0_IK,lmsize) > 0_IK) then
                ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold:ubold-1), getRepeat(rmfill,rmsize)]
            else
                ArrayCentered_ref(lbnew:ubnew) = [Array(lbold:ubold-1)]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array(:) = [.true._LKG, .true._LKG, .true._LKG, .true._LKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold:ubold-1), getRepeat(rmfill,rmsize)]
#elif       IK_ENABLED && D1_ENABLED
            Array(:) = [1_IKG, 2_IKG, 3_IKG, 4_IKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold:ubold-1), getRepeat(rmfill,rmsize)]
#elif       CK_ENABLED && D1_ENABLED
            Array(:) = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (4._CKG,-4._CKG)]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold:ubold-1), getRepeat(rmfill,rmsize)]
#elif       RK_ENABLED && D1_ENABLED
            Array(:) = [1._RKG, 2._RKG, 3._RKG, 4._RKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold:ubold-1), getRepeat(rmfill,rmsize)]
#endif
            cbeg = lbold + getOption(0_IK,lmsize) + abs(sizecentered - sizeold) / 2_IK
            cend = cbeg + min(sizeold,sizecentered) - 1_IK

#if         setCentered_ENABLED
            if (menabled) then
                call setCentered(ArrayCentered, Array, lmsize, rmsize, fill, lmfill, rmfill)
            else
                call setCentered(ArrayCentered, Array, fill)
            end if

            assertion = assertion .and. GET_LBOUND(ArrayCentered) == lbnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the lower bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. GEN_UBOUND(ArrayCentered) == ubnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))
#elif       getCentered_ENABLED
            if (menabled) then
                ArrayCentered = getCentered(Array, sizecentered, lmsize, rmsize, fill, lmfill, rmfill)
            else
                ArrayCentered = getCentered(Array, sizecentered, fill)
            end if
#else
#error      "Unrecognized interface."
#endif

            assertion = assertion .and. GET_SIZE(ArrayCentered) == sizenew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must yield an array of proper size, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            if (present(fill)) then
            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly fill the new elements with `fill` = "//getStr(present(fill)), int(__LINE__, IK))
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Shrink and center array perfectly
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

            lbold = GEN_LBOLD(1_IK)
            ubold = GEN_UBOLD(4_IK)
            lbnew = GEN_LBNEW(1_IK)
            ubnew = GEN_UBNEW(2_IK)
            sizeold = ubold - lbold + 1_IK
            sizenew = ubnew - lbnew + 1_IK
            sizecentered = sizenew - getOption(0_IK,lmsize) - getOption(0_IK,rmsize)
            ALLOC_ARRAY(Array,lbold,ubold)
            ALLOC_ARRAY(ArrayCentered,lbnew,ubnew)
            ALLOC_ARRAY(ArrayCentered_ref,lbnew,ubnew)
#if         SK_ENABLED && D0_ENABLED
            call setUnifRand(Array, repeat(SKG_"A",len(Array)), repeat(SKG_"Z",len(Array)))
            ArrayCentered_ref(lbnew:ubnew) = getRepeat(lmfill,lmsize)//Array(lbold+1:ubold-1)//getRepeat(rmfill,rmsize)
#elif       SK_ENABLED && D1_ENABLED
            Array(:) = [SKG_"AA", SKG_"BB", SKG_"CC", SKG_"DD"]
            if (getOption(0_IK,lmsize) > 0_IK) then
                ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold+1:ubold-1), getRepeat(rmfill,rmsize)]
            else
                ArrayCentered_ref(lbnew:ubnew) = [Array(lbold+1:ubold-1)]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array(:) = [.true._LKG, .true._LKG, .true._LKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold+1:ubold-1), getRepeat(rmfill,rmsize)]
#elif       IK_ENABLED && D1_ENABLED
            Array(:) = [1_IKG, 2_IKG, 3_IKG, 4_IKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold+1:ubold-1), getRepeat(rmfill,rmsize)]
#elif       CK_ENABLED && D1_ENABLED
            Array(:) = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (4._CKG,-4._CKG)]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold+1:ubold-1), getRepeat(rmfill,rmsize)]
#elif       RK_ENABLED && D1_ENABLED
            Array(:) = [1._RKG, 2._RKG, 3._RKG, 4._RKG]
            ArrayCentered_ref(lbnew:ubnew) = [getRepeat(lmfill,lmsize), Array(lbold+1:ubold-1), getRepeat(rmfill,rmsize)]
#endif
            cbeg = lbold + abs(sizecentered - sizeold) / 2_IK
            cend = cbeg + min(sizeold,sizecentered) - 1_IK

#if         setCentered_ENABLED
            if (menabled) then
                call setCentered(ArrayCentered, Array, lmsize, rmsize, fill, lmfill, rmfill)
            else
                call setCentered(ArrayCentered, Array, fill)
            end if


            assertion = assertion .and. GET_LBOUND(ArrayCentered) == lbnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the lower bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. GEN_UBOUND(ArrayCentered) == ubnew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))
#elif       getCentered_ENABLED
            if (menabled) then
                ArrayCentered = getCentered(Array, sizecentered, lmsize, rmsize, fill, lmfill, rmfill)
            else
                ArrayCentered = getCentered(Array, sizecentered, fill)
            end if
#else
#error      "Unrecognized interface."
#endif

            assertion = assertion .and. GET_SIZE(ArrayCentered) == sizenew
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must yield an array of proper size, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly set the upper bound of the output array, with present(fill) = "//getStr(present(fill)), int(__LINE__, IK))

            if (present(fill)) then
            assertion = assertion .and. ALL(ArrayCentered(GET_LBOUND(ArrayCentered)+getOption(0_IK,lmsize):GEN_UBOUND(ArrayCentered)-getOption(0_IK,rmsize)) IS_EQUAL ArrayCentered_ref(lbnew+getOption(0_IK,lmsize):ubnew-getOption(0_IK,rmsize)))
            call report()
            call test%assert(assertion, SK_"Call to setCentered()/getCentered() must properly fill the new elements with `fill` = "//getStr(present(fill)), int(__LINE__, IK))
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
            if (allocated(ArrayCentered)) deallocate(ArrayCentered)
            if (allocated(ArrayCentered_ref)) deallocate(ArrayCentered_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
        PURE function getRepeat(fill, count) result(Array)
            integer(IK)     , intent(in), optional  :: count
            character(1,SKG), intent(in), optional  :: fill
            character(:,SKG), allocatable           :: Array
            Array = repeat(getOption(SKG_" ",fill), getOption(0_IK,count))
        end function
#else
        PURE function getRepeat(fill, count) result(Array)
            integer(IK)     , intent(in), optional  :: count
#if         SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in), optional  :: fill
            character(2,SKG), allocatable           :: Array(:)
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in), optional  :: fill
            logical(LKG)    , allocatable           :: Array(:)
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in), optional  :: fill
            integer(IKG)    , allocatable           :: Array(:)
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in), optional  :: fill
            complex(CKG)    , allocatable           :: Array(:)
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in), optional  :: fill
            real(RKG)       , allocatable           :: Array(:)
#else
#error      "Unrecognized interface."
#endif
            integer(IK) :: i
            Array = [( getOption(fill_def,fill), i = 1_IK, getOption(0_IK,count) )]
        end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array          ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayCentered  ", ArrayCentered
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALLOC_ARRAY_CENTERED
#undef ALLOC_ARRAY
#undef GEN_UBOUND
#undef GET_LBOUND
#undef GEN_UBOLD
#undef GEN_UBNEW
#undef GEN_LBOLD
#undef GEN_LBNEW
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL